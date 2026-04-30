module FixedJ

# v7 §20 — bisect V vs SHE so the model's total Faradaic current density
# matches a target j (mA/cm² → A/m² conversion happens at the call site).
# v7 change: 4-reaction kinetics (added TCH); FixedJResult carries j4/FE_TCH_pct.

using ..Params
using ..Mesh
using ..Chemistry
using ..Kinetics
using ..Assembly
using ..Solver

export solve_at_j, FixedJResult

# ---------- Galvanostatic helper (concrete wrapper for the fixed-j BVP) ----------
function _solve_galvanostatic!(u::Vector{Float64}, V0::Float64, j_target::Float64,
                               mesh, eps_org::Float64, c_eq;
                               max_iter::Int = 60,
                               tol_rel::Float64 = 1.0e-3,
                               verbose::Bool = false)
    build_res = (Vv) -> (F, x) -> full_residual!(F, x, mesh, eps_org, Vv, 1.0, 1.0, c_eq)
    j_fn      = (x, Vv) -> first(_j_total(x, Vv))
    res_g = newton_galvanostatic!(u, V0, j_target, build_res, j_fn;
                                  tol_rel_j = tol_rel,
                                  max_iter  = max_iter,
                                  verbose   = verbose)
    if !res_g.converged
        return res_g
    end
    # Galvanostatic Newton matched j but physics may not be fully converged.
    # Run V-fixed Newton at the galvanostatic V to converge physics from the warm start.
    V_star = res_g.V
    res_v = _solve_at_V!(u, mesh, eps_org, V_star, c_eq;
                          max_iter = max_iter, jacobian_mode = :fd, verbose = false)
    if res_v.converged
        j_final = first(_j_total(u, V_star))
        if abs(j_final - j_target) ≤ tol_rel * j_target
            return (converged = true, iter = res_g.iter, normF = res_v.normF,
                    u = u, V = V_star)
        end
    end
    # V-fixed physics refinement failed or j drifted — return galvanostatic result as-is
    return res_g
end

# ---------- Result type ----------
struct FixedJResult
    converged::Bool
    V_cathode::Float64
    j_total::Float64
    j1::Float64               # ADPN partial current density
    j2::Float64               # PN partial current density
    j3::Float64               # HER partial current density
    j4::Float64               # TCH partial current density
    FE_ADN_pct::Float64
    FE_PN_pct::Float64
    FE_HER_pct::Float64
    FE_TCH_pct::Float64
    state::Vector{Float64}    # converged DOF vector
    n_bisect::Int
    note::String
end

# ---------- Newton solve at one V using a warm start ----------
function _solve_at_V!(u::Vector{Float64}, mesh, eps_org::Float64,
                      V::Float64, c_eq;
                      max_iter::Int = 60,
                      jacobian_mode::Symbol = :ad,
                      verbose::Bool = false)
    residual! = (F, x) -> full_residual!(F, x, mesh, eps_org, V,
                                         1.0, 1.0, c_eq)
    return newton_solve!(u, residual!;
                         max_iter      = max_iter,
                         jacobian_mode = jacobian_mode,
                         verbose       = verbose)
end

# ---------- Total Faradaic current density at the converged state ----------
function _j_total(u::Vector{Float64}, V::Float64)
    j1, j2, j3, j4 = faradaic_currents_from_state(u, V, 1.0)
    return (j1 + j2 + j3 + j4, j1, j2, j3, j4)
end

@inline function _fe(jt, j1, j2, j3, j4)
    FE_ADN = 100.0 * j1 / jt
    FE_PN  = 100.0 * j2 / jt
    FE_HER = 100.0 * j3 / jt
    FE_TCH = 100.0 * j4 / jt
    return FE_ADN, FE_PN, FE_HER, FE_TCH
end

"""
    solve_at_j(j_target_A_m2, eps_org, delta_m, mesh, u_warm, c_eq;
               V_hi=-0.8, V_lo=-2.5, ds_init=0.05, ds_min=5e-4, ds_max=0.20,
               tol_rel=1e-3, max_walk=200, max_bisect=40,
               jacobian_mode=:fd, newton_max_iter=60,
               j0=nothing, alpha_c=nothing, n_orders=nothing,
               V_hint=NaN, u_hint=nothing,
               verbose=false) -> FixedJResult

Find V vs SHE such that the model's total cathodic current density matches
`j_target_A_m2` to within `tol_rel · j_target`.

If `j0` and `alpha_c` are provided as 4-tuples and `n_orders` as a 3-tuple,
they are pushed into the Kinetics override Ref for the duration of this call.

If `V_hint` and `u_hint` are provided (previous solution V_cathode and state),
a fast path is attempted: Newton at V_hint from u_hint. If that converges and
|j - j_target| ≤ tol_rel·j_target, return immediately without any walk. If j
doesn't match (parameter changed more than a tiny FD step), the walk starts from
(V_hint + 0.05) using the converged state at V_hint as warm start.
"""
function solve_at_j(j_target_A_m2::Float64,
                    eps_org::Float64,
                    _delta_m::Float64,
                    mesh,
                    u_warm::Vector{Float64},
                    c_eq;
                    V_hi::Float64       = -0.8,
                    V_lo::Float64       = -2.5,
                    ds_init::Float64    = 0.05,
                    ds_min::Float64     = 5.0e-4,
                    ds_max::Float64     = 0.20,
                    tol_rel::Float64    = 1.0e-3,
                    max_walk::Int       = 200,
                    max_bisect::Int     = 40,
                    jacobian_mode::Symbol = :fd,
                    newton_max_iter::Int  = 60,
                    j0::Union{Nothing,NTuple{4,Float64}}       = nothing,
                    alpha_c::Union{Nothing,NTuple{4,Float64}}  = nothing,
                    n_orders::Union{Nothing,NTuple{3,Float64}} = nothing,
                    V_hint::Float64                            = NaN,
                    u_hint::Union{Nothing,Vector{Float64}}     = nothing,
                    verbose::Bool       = false)

    @assert V_lo < V_hi "V_lo must be more negative than V_hi"
    @assert j_target_A_m2 > 0 "j_target must be positive (cathodic)"

    use_override = !(j0 === nothing || alpha_c === nothing)
    prev_override = Kinetics.KIN_OVERRIDE[]
    if use_override
        n_use = n_orders === nothing ? (2.0, 1.0, 3.0) : n_orders
        set_kinetic_override!(j0, alpha_c, n_use)
    end

    try
        # ---------- Fast path: direct Newton at V_hint ----------
        # When V_hint and u_hint are provided (previous converged V_cathode + state),
        # try Newton at V_hint first. For tiny FD perturbations the solution barely
        # moves, so Newton converges in 1-2 iterations and j ≈ j_target immediately.
        u_walk_warm = u_warm      # warm start used below for the main walk
        V_hi_eff    = V_hi        # V_hi used for the main walk
        if !isnan(V_hint) && u_hint !== nothing
            u_try = copy(u_hint)
            res_h = _solve_at_V!(u_try, mesh, eps_org, V_hint, c_eq;
                                 max_iter = newton_max_iter,
                                 jacobian_mode = jacobian_mode)
            if res_h.converged
                jt, j1, j2, j3, j4 = _j_total(u_try, V_hint)
                if abs(jt - j_target_A_m2) ≤ tol_rel * j_target_A_m2
                    fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1, j2, j3, j4)
                    return FixedJResult(true, V_hint, jt, j1, j2, j3, j4,
                                        fe_adn, fe_pn, fe_her, fe_tch,
                                        u_try, 0, "hint")
                end
                # j doesn't match (theta changed more than an FD step).
                # Start the walk from just above V_hint with u_try as warm start.
                V_hi_eff    = min(V_hint + 0.10, V_hi)
                u_walk_warm = u_try
            end
        end

        # ---------- Step 1: warm-start at V_hi_eff ----------
        u_curr = copy(u_walk_warm)
        res_hi = _solve_at_V!(u_curr, mesh, eps_org, V_hi_eff, c_eq;
                              max_iter = newton_max_iter,
                              jacobian_mode = jacobian_mode,
                              verbose = verbose)
        if !res_hi.converged
            # Warm-start Newton failed at V_hi. This can happen when the cached
            # warm-start state was written during a Jacobian FD perturbation and
            # is inconsistent with the current kinetics. Retry with a cold
            # equilibrium guess at V_hi=-0.3 (near-equilibrium, always converges).
            u_curr = make_initial_guess(length(u_walk_warm) ÷ 10, c_eq, eps_org)
            V_hi_eff = -0.3
            res_hi = _solve_at_V!(u_curr, mesh, eps_org, V_hi_eff, c_eq;
                                  max_iter = newton_max_iter,
                                  jacobian_mode = jacobian_mode,
                                  verbose = verbose)
            if !res_hi.converged
                return FixedJResult(false, V_hi_eff, NaN, NaN, NaN, NaN, NaN,
                                    NaN, NaN, NaN, NaN, u_warm, 0,
                                    "Newton failed at V_hi=$(V_hi_eff) even with cold start")
            end
        end
        j_curr, _, _, _, _ = _j_total(u_curr, V_hi_eff)
        V_curr = V_hi_eff

        if j_curr ≥ j_target_A_m2
            return FixedJResult(false, V_hi_eff, j_curr, NaN, NaN, NaN, NaN,
                                NaN, NaN, NaN, NaN, u_curr, 0,
                                "j(V_hi)=$(j_curr) already ≥ j_target=$(j_target_A_m2); raise V_hi")
        end

        # ---------- Step 2: walk V_curr → V_lo until j ≥ j_target ----------
        ds = ds_init
        u_prev  = copy(u_curr)
        V_prev  = V_curr
        crossed = false
        n_steps = 0
        has_prev = false   # true after the first successful V step

        for step in 1:max_walk
            n_steps = step
            V_next = max(V_curr - ds, V_lo)

            # First-order Euler predictor: extrapolate u along the solution curve.
            # Provides a much better Newton initial guess in the stiff phosphate-buffer
            # regime near V≈-1.63 where the plain warm-start (copy of u_curr) stalls.
            if has_prev && abs(V_curr - V_prev) > 1e-12
                slope   = (V_next - V_curr) / (V_curr - V_prev)
                u_trial = u_curr .+ slope .* (u_curr .- u_prev)
            else
                u_trial = copy(u_curr)
            end
            res     = _solve_at_V!(u_trial, mesh, eps_org, V_next, c_eq;
                                   max_iter = newton_max_iter,
                                   jacobian_mode = jacobian_mode,
                                   verbose = verbose)

            if !res.converged
                ds *= 0.3
                if ds < ds_min
                    # V-walk stalled — fall back to galvanostatic Newton starting from
                    # the last converged state (u_curr, V_curr, j_curr < j_target).
                    verbose && @warn "V-walk stalled; trying galvanostatic fallback" V=V_curr j=j_curr j_target=j_target_A_m2
                    u_gal = copy(u_curr)
                    res_g = _solve_galvanostatic!(u_gal, V_curr, j_target_A_m2,
                                                  mesh, eps_org, c_eq;
                                                  max_iter = newton_max_iter,
                                                  tol_rel  = tol_rel,
                                                  verbose  = verbose)
                    if res_g.converged
                        jt, j1, j2, j3, j4 = _j_total(res_g.u, res_g.V)
                        fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1, j2, j3, j4)
                        return FixedJResult(true, res_g.V, jt, j1, j2, j3, j4,
                                            fe_adn, fe_pn, fe_her, fe_tch,
                                            res_g.u, step, "galvanostatic fallback")
                    end

                    # Second attempt: start from V_lo with AN-depleted initial guess.
                    # When j_target > AN transport limit, the physical solution has
                    # c_AN→0 at the electrode with HER carrying most of the current.
                    # This basin is unreachable from the AN-rich V_hi end but accessible
                    # by starting at V_lo where HER already dominates.
                    verbose && @warn "galvanostatic from AN-rich end failed; trying from V_lo with depleted AN guess"
                    N_dep  = length(u_warm) ÷ 10
                    u_dep  = copy(u_warm)
                    c_bulk = Chemistry.c_AN_bulk(eps_org)
                    c_seed = 1.0e-3
                    for ix in 1:N_dep
                        frac = (ix - 1) / max(1, N_dep - 1)   # 0=electrode, 1=outer
                        c_AN_ix = c_seed + frac * (c_bulk - c_seed)
                        u_dep[10*(ix-1) + 6] = log(max(c_AN_ix, c_seed))
                    end
                    res_dep = _solve_at_V!(u_dep, mesh, eps_org, V_lo, c_eq;
                                           max_iter = newton_max_iter,
                                           jacobian_mode = jacobian_mode)
                    if res_dep.converged
                        jt_lo, _, _, _, _ = _j_total(u_dep, V_lo)
                        if jt_lo >= j_target_A_m2
                            res_g2 = _solve_galvanostatic!(u_dep, V_lo, j_target_A_m2,
                                                            mesh, eps_org, c_eq;
                                                            max_iter = newton_max_iter,
                                                            tol_rel  = tol_rel,
                                                            verbose  = verbose)
                            if res_g2.converged
                                jt, j1, j2, j3, j4 = _j_total(res_g2.u, res_g2.V)
                                fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1, j2, j3, j4)
                                return FixedJResult(true, res_g2.V, jt, j1, j2, j3, j4,
                                                    fe_adn, fe_pn, fe_her, fe_tch,
                                                    res_g2.u, step,
                                                    "galvanostatic from V_lo (AN-depleted)")
                            end
                        end
                    end
                    return FixedJResult(false, V_curr, NaN, NaN, NaN, NaN, NaN,
                                        NaN, NaN, NaN, NaN, u_curr, step,
                                        "Walk shrunk below ds_min=$(ds_min) at V≈$(V_curr); " *
                                        "both galvanostatic fallbacks failed")
                end
                continue
            end

            j_next, _, _, _, _ = _j_total(u_trial, V_next)

            u_prev = u_curr; V_prev = V_curr
            u_curr = u_trial; V_curr = V_next; j_curr = j_next
            has_prev = true

            if j_next ≥ j_target_A_m2
                crossed = true
                break
            end

            if V_curr ≤ V_lo + 1e-12
                return FixedJResult(false, V_lo, j_curr, NaN, NaN, NaN, NaN,
                                    NaN, NaN, NaN, NaN, u_curr, step,
                                    "Reached V_lo=$(V_lo) with j=$(j_curr) < j_target=$(j_target_A_m2)")
            end

            if res.iter ≤ 4
                ds = min(ds * 1.4, ds_max)
            elseif res.iter ≤ 10
                ds = min(ds * 1.1, ds_max)
            elseif res.iter > 20
                ds = max(ds * 0.7, ds_min)
            end
        end

        if !crossed
            return FixedJResult(false, V_curr, j_curr, NaN, NaN, NaN, NaN,
                                NaN, NaN, NaN, NaN, u_curr, n_steps,
                                "max_walk=$(max_walk) exhausted before crossing j_target")
        end

        # Bracket: a (more negative, j ≥ target) < b (less negative, j < target)
        a, b = V_curr, V_prev
        u_a, u_b = u_curr, u_prev

        # Early-out if already within tolerance
        if abs(j_curr - j_target_A_m2) ≤ tol_rel * j_target_A_m2
            jt, j1, j2, j3, j4 = _j_total(u_curr, V_curr)
            fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1, j2, j3, j4)
            return FixedJResult(true, V_curr, jt, j1, j2, j3, j4,
                                fe_adn, fe_pn, fe_her, fe_tch, u_curr, 0, "")
        end

        # ---------- Step 3: bisection within (a, b) ----------
        u_mid = copy(u_a)
        V_mid = 0.5 * (a + b)
        jt = j1m = j2m = j3m = j4m = NaN

        for k in 1:max_bisect
            V_mid  = 0.5 * (a + b)
            u_seed = abs(V_mid - a) < abs(V_mid - b) ? u_a : u_b
            u_mid  = copy(u_seed)
            res_m  = _solve_at_V!(u_mid, mesh, eps_org, V_mid, c_eq;
                                  max_iter = newton_max_iter,
                                  jacobian_mode = jacobian_mode,
                                  verbose = verbose)
            if !res_m.converged
                return FixedJResult(false, V_mid, NaN, NaN, NaN, NaN, NaN,
                                    NaN, NaN, NaN, NaN, u_a, k,
                                    "Newton failed inside bisection at V_mid=$(V_mid) (step $k)")
            end
            jt, j1m, j2m, j3m, j4m = _j_total(u_mid, V_mid)
            if abs(jt - j_target_A_m2) ≤ tol_rel * j_target_A_m2
                fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1m, j2m, j3m, j4m)
                return FixedJResult(true, V_mid, jt, j1m, j2m, j3m, j4m,
                                    fe_adn, fe_pn, fe_her, fe_tch, u_mid, k, "")
            end
            if jt > j_target_A_m2
                a = V_mid; u_a = u_mid
            else
                b = V_mid; u_b = u_mid
            end
        end

        fe_adn, fe_pn, fe_her, fe_tch = _fe(jt, j1m, j2m, j3m, j4m)
        return FixedJResult(false, V_mid, jt, j1m, j2m, j3m, j4m,
                            fe_adn, fe_pn, fe_her, fe_tch, u_mid, max_bisect,
                            "Bisection did not reach tol_rel=$(tol_rel); " *
                            "|Δj/j| = $(abs(jt - j_target_A_m2)/j_target_A_m2)")
    finally
        if use_override
            Kinetics.KIN_OVERRIDE[] = prev_override
        end
    end
end

end # module
