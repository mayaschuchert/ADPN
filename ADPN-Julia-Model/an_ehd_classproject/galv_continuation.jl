module GalvContinuation

# v7 §21 — Galvanostatic continuation solver. Replaces fixed_j_solver's
# V_hi-heuristic V-walk + bisection with anchor-and-step in log10(j_total).
#
# Adapted from SOLVER_STACK.md (CO2R Phase-15) with our system's specifics:
#  - Anchor at small j_anchor (≈ 1 A/m²) where Newton trivially converges
#    from the equilibrium guess at V ≈ -0.5 V (η_HER = +0.33 V → no kinetics).
#  - Continuation parameter: log10(j_total).  Uniform log-j stepping
#    auto-adapts to the Tafel exponential (every decade ≈ RT/(αF) volts).
#  - Augmented (n+1)×(n+1) Newton at each step jointly updates (u, V) so
#    j_total(u, V) hits the target — no V-bisection needed.
#  - Warm-start fast path: if a (u, V) hint is provided AND |Δlog10(j)| is
#    small, try a single galvanostatic Newton at the hint and skip the chain.
#  - Fallback: full chain from the anchor.  Robust across kinetics changes.
#
# This module never reads from / writes to the FitContext warm-start cache;
# the caller decides whether to update its cache after a successful return.
# That isolation prevents the Jacobian-FD warm-start corruption that broke
# fixed_j_solver under LM Jacobian evaluations.

using Printf
using ..Params
using ..Mesh
using ..Chemistry
using ..Kinetics
using ..Assembly
using ..Solver

export ContResult, solve_at_j_cont

# ------------------------------------------------------------
# Result type — mirrors FixedJResult so callers (run_stage4*) can swap.
# ------------------------------------------------------------
struct ContResult
    converged::Bool
    V_cathode::Float64
    j_total::Float64
    j1::Float64               # ADPN partial current density [A m⁻²]
    j2::Float64               # PN partial current density
    j3::Float64               # HER partial current density
    j4::Float64               # TCH partial current density
    FE_ADN_pct::Float64
    FE_PN_pct::Float64
    FE_HER_pct::Float64
    FE_TCH_pct::Float64
    state::Vector{Float64}
    n_newton::Int             # total Newton iterations across all log-j steps
    note::String
end

# ------------------------------------------------------------
# Helpers (private)
# ------------------------------------------------------------
@inline function _j_total(u::Vector{Float64}, V::Float64)
    j1, j2, j3, j4 = faradaic_currents_from_state(u, V, 1.0)
    return (j1 + j2 + j3 + j4, j1, j2, j3, j4)
end

@inline function _fe(jt, j1, j2, j3, j4)
    100.0 * j1 / jt, 100.0 * j2 / jt, 100.0 * j3 / jt, 100.0 * j4 / jt
end

# Newton at fixed V (used to (a) converge the anchor and (b) compute j(V)).
function _solve_at_V!(u::Vector{Float64}, mesh, eps_org::Float64,
                      V::Float64, c_eq;
                      max_iter::Int = 60,
                      jacobian_mode::Symbol = :fd,
                      verbose::Bool = false)
    residual! = (F, x) -> full_residual!(F, x, mesh, eps_org, V,
                                         1.0, 1.0, c_eq)
    return newton_solve!(u, residual!;
                         max_iter      = max_iter,
                         jacobian_mode = jacobian_mode,
                         verbose       = verbose)
end

# Bordered (n+1)×(n+1) galvanostatic Newton at a single j target.  Wraps
# Solver.newton_galvanostatic! with our residual builder and j_total fn.
function _galv_newton!(u::Vector{Float64}, V0::Float64, j_target::Float64,
                       mesh, eps_org::Float64, c_eq;
                       max_iter::Int = 40,
                       tol_rel_j::Float64 = 1.0e-3,
                       verbose::Bool = false)
    build_res = (Vv) -> (F, x) -> full_residual!(F, x, mesh, eps_org, Vv,
                                                 1.0, 1.0, c_eq)
    j_fn      = (x, Vv) -> first(_j_total(x, Vv))
    return newton_galvanostatic!(u, V0, j_target, build_res, j_fn;
                                 tol_rel_j = tol_rel_j,
                                 max_iter  = max_iter,
                                 verbose   = verbose)
end

# Build the equilibrium anchor: cold equilibrium guess + Newton at V_anchor.
# Stepping V_anchor more cathodic in 50 mV increments until j ≥ j_floor.
# Returns (u_anchor, V_anchor, j_anchor, n_iter, ok).
function _build_anchor(mesh, eps_org::Float64, c_eq;
                       V_anchor_init::Float64 = -0.5,
                       j_floor::Float64       = 0.1,    # A/m² (= 0.01 mA/cm²)
                       V_floor::Float64       = -1.6,   # past HER turn-on but pre-ADN/PN
                       dV::Float64            = -0.05,
                       max_iter::Int          = 60)
    N = length(mesh.dx)
    u = make_initial_guess(N, c_eq, eps_org)
    V = V_anchor_init
    n_iter_total = 0

    res = _solve_at_V!(u, mesh, eps_org, V, c_eq; max_iter = max_iter)
    n_iter_total += res.iter
    if !res.converged
        return (u = u, V = V, j = NaN, n = n_iter_total, ok = false,
                msg = "anchor Newton failed at V=$V_anchor_init")
    end

    j, _, _, _, _ = _j_total(u, V)
    # Ramp V more cathodic until j > j_floor (so log10 is well defined).
    while j < j_floor && V > V_floor
        V_next = max(V + dV, V_floor)
        res = _solve_at_V!(u, mesh, eps_org, V_next, c_eq; max_iter = max_iter)
        n_iter_total += res.iter
        if !res.converged
            return (u = u, V = V, j = j, n = n_iter_total, ok = false,
                    msg = "anchor ramp Newton failed at V=$V_next")
        end
        V = V_next
        j, _, _, _, _ = _j_total(u, V)
    end

    if j < j_floor
        return (u = u, V = V, j = j, n = n_iter_total, ok = false,
                msg = "anchor ramp reached V_floor=$V_floor with j=$j < j_floor=$j_floor")
    end
    return (u = u, V = V, j = j, n = n_iter_total, ok = true, msg = "")
end

# ------------------------------------------------------------
# Main entry: galvanostatic continuation from anchor → j_target
# ------------------------------------------------------------
"""
    solve_at_j_cont(j_target, eps_org, mesh, c_eq;
                    j0=nothing, alpha_c=nothing, n_orders=nothing,
                    fixed_j=nothing,
                    warm_u=nothing, warm_V=NaN, warm_j=NaN,
                    n_per_decade=8, fast_logj_window=0.05,
                    tol_rel_j=1e-3, newton_max_iter=40, verbose=false)
        -> ContResult

Galvanostatic continuation in log10(j_total).

Mandatory:
- `j_target`   target total current density [A m⁻²]
- `eps_org`    organic-phase volume fraction
- `mesh`, `c_eq`  per-row mesh and bulk equilibrium

Kinetic overrides (optional, for kinetics-fitting use):
- `j0`, `alpha_c` 4-tuples (ADN, PN, HER, TCH) — pushed into Kinetics.KIN_OVERRIDE
- `n_orders`      3-tuple (n_ADN, n_PN, n_TCH)
- `fixed_j`       4-tuple of A/m² (NaN sentinel = unset).  Pushed into Kinetics.FIXED_J
                  — used by Lees S1.11 sequential fitting to clamp branches.

Warm start (optional fast path):
- `warm_u`, `warm_V`, `warm_j` → if |log10(j_target/warm_j)| ≤ fast_logj_window,
  attempt a single galvanostatic Newton at (warm_u, warm_V).  On failure, fall
  back to the full anchor → log-j chain.
"""
function solve_at_j_cont(j_target_A_m2::Float64,
                         eps_org::Float64,
                         mesh,
                         c_eq;
                         j0::Union{Nothing,NTuple{4,Float64}}       = nothing,
                         alpha_c::Union{Nothing,NTuple{4,Float64}}  = nothing,
                         n_orders::Union{Nothing,NTuple{3,Float64}} = nothing,
                         fixed_j::Union{Nothing,NTuple{4,Float64}}  = nothing,
                         warm_u::Union{Nothing,Vector{Float64}}     = nothing,
                         warm_V::Float64                            = NaN,
                         warm_j::Float64                            = NaN,
                         n_per_decade::Int                          = 8,
                         fast_logj_window::Float64                  = 0.05,
                         tol_rel_j::Float64                         = 1.0e-3,
                         newton_max_iter::Int                       = 40,
                         verbose::Bool                              = false)

    @assert j_target_A_m2 > 0 "j_target must be positive (cathodic)"

    # ---- push optional kinetic overrides for the duration of this call ----
    # KIN_OVERRIDE (j0/αc/n) applies throughout. FIXED_J (clamping for S1.11)
    # is held off during anchor build — clamping a branch high makes anchoring
    # at j≈0 impossible — and re-enabled before the log-j continuation chain.
    use_kin_override = !(j0 === nothing || alpha_c === nothing)
    prev_kin_override = Kinetics.KIN_OVERRIDE[]
    prev_fixed_j      = Kinetics.FIXED_J[]
    if use_kin_override
        n_use = n_orders === nothing ? (2.0, 1.0, 3.0) : n_orders
        set_kinetic_override!(j0, alpha_c, n_use)
    end
    # NOTE: FIXED_J intentionally NOT set yet — defer until after anchor.

    n_iter_total = 0
    try
        # ---- Fast path: single galvanostatic Newton from a warm hint ----
        # Warm hint is at the same j_target → clamping is consistent. Enable
        # FIXED_J and try one Newton.
        if warm_u !== nothing && !isnan(warm_V) && !isnan(warm_j) && warm_j > 0.0
            dlogj = abs(log10(j_target_A_m2) - log10(warm_j))
            if dlogj ≤ fast_logj_window
                fixed_j !== nothing && (Kinetics.FIXED_J[] = fixed_j)
                u_try = copy(warm_u)
                res = _galv_newton!(u_try, warm_V, j_target_A_m2,
                                    mesh, eps_org, c_eq;
                                    max_iter = newton_max_iter,
                                    tol_rel_j = tol_rel_j)
                n_iter_total += res.iter
                if res.converged
                    jt, j1, j2, j3, j4 = _j_total(u_try, res.V)
                    fa, fp, fh, ft = _fe(jt, j1, j2, j3, j4)
                    return ContResult(true, res.V, jt, j1, j2, j3, j4,
                                      fa, fp, fh, ft, u_try, n_iter_total,
                                      "fast (Δlog10j=$(round(dlogj; digits=3)))")
                end
                # Roll back FIXED_J for the anchor build below
                Kinetics.FIXED_J[] = prev_fixed_j
            end
        end

        # ---- Build anchor (NO clamping — needs to find j≈0 state) ----
        anc = _build_anchor(mesh, eps_org, c_eq; max_iter = newton_max_iter)
        n_iter_total += anc.n
        if !anc.ok
            return ContResult(false, NaN, NaN, NaN, NaN, NaN, NaN,
                              NaN, NaN, NaN, NaN, Float64[], n_iter_total,
                              "anchor failed: $(anc.msg)")
        end
        u  = anc.u
        V  = anc.V
        j  = anc.j

        # ---- Strategy split: clamped vs unclamped ----
        # When clamping is active, the FIXED_J floor j_clamp_sum can be near
        # j_target (e.g., S1.11 fitting HER while ADN/PN/TCH clamps account
        # for ~88% of j_total).  In that regime the previous "ramp to a
        # floor with headroom, then chain in clamped mode" fails because the
        # un-clamped ramp can't smoothly reach a state where the clamps make
        # sense.  Instead:
        #
        #   (1) Run the FULL un-clamped log-j chain from anchor → j_target,
        #       producing a converged (u, V) at the correct total current
        #       under the un-clamped kinetics.
        #   (2) Enable FIXED_J, then do ONE galvanostatic Newton at j_target
        #       starting from the un-clamped solution.  The Newton needs only
        #       to redistribute partial currents; j_total is already correct.
        #
        # When NOT clamping, behaviour is unchanged — the chain runs once.
        log10_j_target = log10(j_target_A_m2)
        log10_j_anchor = log10(j)
        n_decades = abs(log10_j_target - log10_j_anchor)
        N_steps   = max(2, Int(round(n_decades * n_per_decade)))
        dlogj     = (log10_j_target - log10_j_anchor) / N_steps

        for step in 1:N_steps
            log10_j_step = log10_j_anchor + step * dlogj
            j_step       = 10.0^log10_j_step
            res = _galv_newton!(u, V, j_step, mesh, eps_org, c_eq;
                                max_iter = newton_max_iter,
                                tol_rel_j = tol_rel_j,
                                verbose = false)
            n_iter_total += res.iter
            if !res.converged
                # Try refining: subdivide the failing step
                ok_sub = false
                for sub_div in (2, 4, 8)
                    log10_j_prev = log10_j_anchor + (step - 1) * dlogj
                    sub_dlogj    = dlogj / sub_div
                    u_sub = copy(u);  V_sub = V
                    sub_ok = true
                    for s in 1:sub_div
                        j_sub = 10.0^(log10_j_prev + s * sub_dlogj)
                        rs = _galv_newton!(u_sub, V_sub, j_sub, mesh, eps_org, c_eq;
                                            max_iter = newton_max_iter,
                                            tol_rel_j = tol_rel_j)
                        n_iter_total += rs.iter
                        if !rs.converged
                            sub_ok = false; break
                        end
                        u_sub .= u_sub
                        V_sub  = rs.V
                    end
                    if sub_ok
                        u .= u_sub; V = V_sub
                        ok_sub = true; break
                    end
                end
                if !ok_sub
                    return ContResult(false, V, NaN, NaN, NaN, NaN, NaN,
                                      NaN, NaN, NaN, NaN, u, n_iter_total,
                                      "log-j Newton failed at step $step / $N_steps " *
                                      "(j_step=$(j_step), even with subdivision)")
                end
            else
                V = res.V
            end
        end

        # ---- Phase 2: enable FIXED_J and redistribute partial currents ----
        # Un-clamped chain has produced (u, V) at j_total ≈ j_target. Now
        # clamp the prescribed branches and let one galvanostatic Newton
        # redistribute. Newton step clamps in φ_l (≤ 15 mV) and log-c (≤ 5)
        # bound the redistribution magnitude per iter, so this typically
        # converges in 5–15 iters from the un-clamped warm start.
        if fixed_j !== nothing
            Kinetics.FIXED_J[] = fixed_j
            res_clamp = _galv_newton!(u, V, j_target_A_m2, mesh, eps_org, c_eq;
                                       max_iter  = newton_max_iter * 2,
                                       tol_rel_j = tol_rel_j)
            n_iter_total += res_clamp.iter
            if !res_clamp.converged
                return ContResult(false, V, NaN, NaN, NaN, NaN, NaN,
                                  NaN, NaN, NaN, NaN, u, n_iter_total,
                                  "FIXED_J redistribution failed after un-clamped chain")
            end
            V = res_clamp.V
        end

        jt, j1, j2, j3, j4 = _j_total(u, V)
        fa, fp, fh, ft = _fe(jt, j1, j2, j3, j4)
        note_str = fixed_j === nothing ?
                    "continuation ($N_steps steps)" :
                    "continuation ($N_steps steps) + FIXED_J redistribution"
        return ContResult(true, V, jt, j1, j2, j3, j4,
                          fa, fp, fh, ft, u, n_iter_total, note_str)
    finally
        if use_kin_override
            Kinetics.KIN_OVERRIDE[] = prev_kin_override
        end
        if fixed_j !== nothing
            Kinetics.FIXED_J[] = prev_fixed_j
        end
    end
end

end # module
