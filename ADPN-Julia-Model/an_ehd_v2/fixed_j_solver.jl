module FixedJ

# v6 §20 — bisect V vs SHE so the model's total Faradaic current density
# matches a target j (mA/cm² → A/m² conversion happens at the call site).
# Inside, this wraps the existing v5 Newton continuation: at each candidate
# V the residual is built and Newton is run from a warm start.

using ..Params
using ..Mesh
using ..Chemistry
using ..Kinetics
using ..Assembly
using ..Solver

export solve_at_j, FixedJResult

# ---------- Result type ----------
struct FixedJResult
    converged::Bool
    V_cathode::Float64        # V vs SHE
    j_total::Float64          # A m⁻², achieved
    j1::Float64               # ADPN partial current density
    j2::Float64               # PN partial current density
    j3::Float64               # HER partial current density
    FE_ADN_pct::Float64
    FE_PN_pct::Float64
    FE_HER_pct::Float64
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
    j1, j2, j3 = faradaic_currents_from_state(u, V, 1.0)
    return (j1 + j2 + j3, j1, j2, j3)
end

"""
    solve_at_j(j_target_A_m2, eps_org, delta_m, mesh, u_warm, c_eq;
               V_hi=-0.8, V_lo=-2.5, ds_init=0.05, ds_min=5e-4, ds_max=0.20,
               tol_rel=1e-3, max_walk=200, max_bisect=40,
               jacobian_mode=:ad, newton_max_iter=60,
               j0=nothing, alpha_c=nothing, verbose=false) -> FixedJResult

Find V vs SHE in (V_lo, V_hi) such that the model's total cathodic current
density matches `j_target_A_m2` to within `tol_rel · j_target`.

Algorithm (mirrors v5 `newton_continuation` then bisects within the crossing
interval — no cold V-jump):

  1. Converge at V_hi from `u_warm` (low current; cold-friendly).
  2. Walk V down adaptively in steps `ds`:
       - Each step warm-starts Newton from the previously converged state.
       - On convergence, query j_total. If j_total ≥ j_target, we have crossed
         the target → break with bracket (V_prev, V_step) and converged states.
       - On failure, halve ds and retry.
       - Honour ds_max growth on fast successes (mirrors §10.5 logic).
  3. Bisect inside the crossing interval, warm-starting each Newton from the
     closer of the two bracket states.

`u_warm` should be a v5-bootstrapped equilibrium (or a previous fit row's
converged state). A pristine `make_initial_guess` is fine because step 1's
Newton at V_hi=-0.8 V is essentially zero-current and converges in 1 iter
from the cold guess.

If `j0` and `alpha_c` are provided as 3-tuples, they are pushed into the
Kinetics override Ref for the duration of this call (restored on exit).
`n_orders` (a 2-tuple of AN reaction orders for ADPN and PN) is optional —
if omitted, the override carries the historical defaults `(2.0, 1.0)`.
"""
function solve_at_j(j_target_A_m2::Float64,
                    eps_org::Float64,
                    delta_m::Float64,
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
                    jacobian_mode::Symbol = :ad,
                    newton_max_iter::Int  = 60,
                    j0::Union{Nothing,NTuple{3,Float64}}      = nothing,
                    alpha_c::Union{Nothing,NTuple{3,Float64}} = nothing,
                    n_orders::Union{Nothing,NTuple{2,Float64}} = nothing,
                    verbose::Bool       = false)

    @assert V_lo < V_hi "V_lo must be more negative than V_hi"
    @assert j_target_A_m2 > 0 "j_target must be positive (cathodic)"

    use_override = !(j0 === nothing || alpha_c === nothing)
    prev_override = Kinetics.KIN_OVERRIDE[]
    if use_override
        n_use = n_orders === nothing ? (2.0, 1.0) : n_orders
        set_kinetic_override!(j0, alpha_c, n_use)
    end

    try
        # ---------- Step 1: warm-start at V_hi ----------
        u_curr = copy(u_warm)
        res_hi = _solve_at_V!(u_curr, mesh, eps_org, V_hi, c_eq;
                              max_iter = newton_max_iter,
                              jacobian_mode = jacobian_mode,
                              verbose = verbose)
        if !res_hi.converged
            return FixedJResult(false, V_hi, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_warm, 0, "Newton failed at V_hi=$(V_hi)")
        end
        j_curr, _, _, _ = _j_total(u_curr, V_hi)
        V_curr = V_hi

        # If V_hi already overshoots j_target the bracket is degenerate — bail.
        if j_curr ≥ j_target_A_m2
            return FixedJResult(false, V_hi, j_curr, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_curr, 0,
                                "j(V_hi)=$(j_curr) already ≥ j_target=$(j_target_A_m2); raise V_hi")
        end

        # ---------- Step 2: walk V_curr → V_lo until j ≥ j_target ----------
        ds = ds_init
        u_prev   = copy(u_curr)
        V_prev   = V_curr
        j_prev   = j_curr
        crossed  = false
        n_steps  = 0
        n_fail   = 0

        for step in 1:max_walk
            n_steps = step
            V_next = max(V_curr - ds, V_lo)

            u_trial = copy(u_curr)
            res     = _solve_at_V!(u_trial, mesh, eps_org, V_next, c_eq;
                                   max_iter = newton_max_iter,
                                   jacobian_mode = jacobian_mode,
                                   verbose = verbose)

            if !res.converged
                ds *= 0.3
                n_fail += 1
                if ds < ds_min
                    return FixedJResult(false, V_next, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                        u_curr, step,
                                        "Walk shrunk below ds_min=$(ds_min) at V≈$(V_curr); " *
                                        "Newton failed at V_next=$(V_next)")
                end
                continue
            end

            j_next, _, _, _ = _j_total(u_trial, V_next)

            # advance
            u_prev = u_curr; V_prev = V_curr; j_prev = j_curr
            u_curr = u_trial; V_curr = V_next; j_curr = j_next

            if j_next ≥ j_target_A_m2
                crossed = true
                break
            end

            if V_curr ≤ V_lo + 1e-12
                # reached V_lo without crossing; j(V_lo) < j_target
                return FixedJResult(false, V_lo, j_curr, NaN, NaN, NaN, NaN, NaN, NaN,
                                    u_curr, step,
                                    "Reached V_lo=$(V_lo) with j=$(j_curr) < j_target=$(j_target_A_m2)")
            end

            # adaptive step size — mirror solver.jl §10.5
            if res.iter ≤ 4
                ds = min(ds * 1.4, ds_max)
            elseif res.iter ≤ 10
                ds = min(ds * 1.1, ds_max)
            elseif res.iter > 20
                ds = max(ds * 0.7, ds_min)
            end
        end

        if !crossed
            return FixedJResult(false, V_curr, j_curr, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_curr, n_steps,
                                "max_walk=$(max_walk) exhausted before crossing j_target")
        end

        # Bracket: V_prev (less negative, j_prev < j_target) < V_curr (more
        # negative, j_curr ≥ j_target). Note we are sweeping V more negative,
        # but in numeric terms V_curr < V_prev.
        a, b = V_curr, V_prev      # a < b, j(a) ≥ j_target > j(b)
        u_a, u_b = u_curr, u_prev
        j_a, j_b = j_curr, j_prev

        # Early-out: if j_curr is already within tolerance, we're done.
        if abs(j_curr - j_target_A_m2) ≤ tol_rel * j_target_A_m2
            j_total, j1, j2, j3 = _j_total(u_curr, V_curr)
            FE_ADN = 100 * j1 / j_total
            FE_PN  = 100 * j2 / j_total
            FE_HER = 100 * j3 / j_total
            return FixedJResult(true, V_curr, j_total, j1, j2, j3,
                                FE_ADN, FE_PN, FE_HER, u_curr, 0, "")
        end

        # ---------- Step 3: bisection within (a, b) ----------
        u_mid = copy(u_a)
        V_mid = 0.5 * (a + b)
        j_mid = NaN; j1m = j2m = j3m = NaN

        for k in 1:max_bisect
            V_mid = 0.5 * (a + b)
            u_seed = abs(V_mid - a) < abs(V_mid - b) ? u_a : u_b
            u_mid  = copy(u_seed)
            res_m  = _solve_at_V!(u_mid, mesh, eps_org, V_mid, c_eq;
                                  max_iter = newton_max_iter,
                                  jacobian_mode = jacobian_mode,
                                  verbose = verbose)
            if !res_m.converged
                return FixedJResult(false, V_mid, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                    u_a, k,
                                    "Newton failed inside bisection at V_mid=$(V_mid) (step $k)")
            end
            j_mid, j1m, j2m, j3m = _j_total(u_mid, V_mid)
            if abs(j_mid - j_target_A_m2) ≤ tol_rel * j_target_A_m2
                FE_total = j1m + j2m + j3m
                FE_ADN  = 100 * j1m / FE_total
                FE_PN   = 100 * j2m / FE_total
                FE_HER  = 100 * j3m / FE_total
                return FixedJResult(true, V_mid, j_mid, j1m, j2m, j3m,
                                    FE_ADN, FE_PN, FE_HER, u_mid, k, "")
            end
            if j_mid > j_target_A_m2
                a = V_mid; u_a = u_mid; j_a = j_mid
            else
                b = V_mid; u_b = u_mid; j_b = j_mid
            end
        end

        FE_total = j1m + j2m + j3m
        FE_ADN  = 100 * j1m / FE_total
        FE_PN   = 100 * j2m / FE_total
        FE_HER  = 100 * j3m / FE_total
        return FixedJResult(false, V_mid, j_mid, j1m, j2m, j3m,
                            FE_ADN, FE_PN, FE_HER,
                            u_mid, max_bisect,
                            "Bisection did not reach tol_rel=$(tol_rel); " *
                            "|Δj/j| = $(abs(j_mid - j_target_A_m2)/j_target_A_m2)")
    finally
        if use_override
            Kinetics.KIN_OVERRIDE[] = prev_override
        end
    end
end

end # module
