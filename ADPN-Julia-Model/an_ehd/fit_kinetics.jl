module FitKinetics

# v6 §20 — kinetics-only fit on Bloomquist data with all transport frozen.
#
#   theta = (j0_1, j0_2, j0_3, alpha_c1, alpha_c2, alpha_c3)
#
# Row filters (§20.1):
#   Core     : gap ∈ {0.5, 1.0} mm AND j ≤ 190 mA/cm² AND ε_org ≥ 0.04   (≈60 rows)
#   Extended : gap ∈ {0.5, 1.0} mm AND ε_org ≥ 0.04                       (≈96 rows)
#   Holdout  : gap = 0.25 mm        AND ε_org ≥ 0.04                       (≈48 rows)
#
# Stage 4a fits θ on Core. Stage 4b applies θ forward to Extended and Holdout
# without re-fitting. V_CE and R_contact are frozen at v6 defaults (§20.5).
#
# This module supplies:
#   * row-selection filters returning row index vectors into the master CSV
#   * a transport pre-cache (delta_lev, kappa_eff per unique (gap, Q_total, eps_org))
#   * a residual builder F(theta) that returns the per-row (FE_ADN_resid, FE_PN_resid)
#   * a minimal pure-Julia Levenberg-Marquardt driver (no external deps)
#
# Heavy IO (CSV read, plotting) lives in run_stage4.jl, not here.

using Printf
using LinearAlgebra
using ..Params
using ..Mesh
using ..Chemistry
using ..Hydrodynamics
using ..CellVoltage
using ..FixedJ

export BloomquistRow, FitContext, build_context,
       select_core, select_extended, select_holdout,
       residuals!, loss,
       theta_to_physical, physical_to_theta,
       lm_fit, LMResult,
       N_THETA, THETA_LB, THETA_UB, THETA0

# ---------- Parameter vector layout ----------
const N_THETA = 6
# theta is in *log10* space for j0,r and *linear* for α_c,r — keeps the LM
# step well-conditioned across the j0 ∈ [10⁻⁶, 10⁻¹] A m⁻² fit range.
const THETA_LB = [-6.0, -6.0, -8.0, 0.30, 0.30, 0.30]      # log10 j0 lo + αc lo
const THETA_UB = [-1.0, -1.0, -3.0, 0.70, 0.70, 0.50]      # log10 j0 hi + αc hi
const THETA0   = [log10(1.0e-3), log10(1.0e-3), log10(1.0e-5),
                  0.50, 0.50, 0.40]

"Convert internal theta vector → physical (j0_tuple, alpha_c_tuple)."
function theta_to_physical(theta::AbstractVector{Float64})
    j0 = (10.0^theta[1], 10.0^theta[2], 10.0^theta[3])
    ac = (theta[4], theta[5], theta[6])
    return j0, ac
end

"Inverse of theta_to_physical for diagnostics / initial guesses."
function physical_to_theta(j0::NTuple{3,Float64}, ac::NTuple{3,Float64})
    return [log10(j0[1]), log10(j0[2]), log10(j0[3]), ac[1], ac[2], ac[3]]
end

# ---------- Row container ----------
"A single Bloomquist row with model-side derived quantities cached."
struct BloomquistRow
    table::String
    gap_mm::Float64
    Q_total_mL_min::Float64
    j_mA_cm2::Float64
    phi_AN::Float64           # ε_org
    Q_aq_mL_min::Float64
    Q_org_mL_min::Float64
    We_aq::Float64
    We_org::Float64
    FE_ADN_pct::Float64       # observed
    FE_TCH_pct::Float64
    FE_PN_pct::Float64
    PR_ADN_kg_cm2_h::Float64
    EP_ADN_kg_kWh::Float64
    # Derived (set by build_context)
    gap_m::Float64
    j_target_A_m2::Float64
    delta_lev_m::Float64      # Lévêque BL thickness
    R_series_Ohm_m2::Float64  # cached for V_cell prediction
end

# ---------- Filters (§20.1) ----------
function select_core(rows::Vector{BloomquistRow})
    return [i for (i, r) in pairs(rows)
            if (r.gap_mm == 0.5 || r.gap_mm == 1.0) &&
               r.j_mA_cm2 ≤ 190.0 &&
               r.phi_AN ≥ 0.04]
end

function select_extended(rows::Vector{BloomquistRow})
    return [i for (i, r) in pairs(rows)
            if (r.gap_mm == 0.5 || r.gap_mm == 1.0) &&
               r.phi_AN ≥ 0.04]
end

function select_holdout(rows::Vector{BloomquistRow})
    return [i for (i, r) in pairs(rows)
            if r.gap_mm == 0.25 &&
               r.phi_AN ≥ 0.04]
end

# ---------- Fit context ----------
"""
Bundle of read-only objects shared across fit-loss evaluations:
- the parsed Bloomquist row vector
- the row-index selection (Core / Extended / Holdout)
- mesh template (one mesh per unique δ; built lazily)
- bulk equilibrium (shared — composition is the same for all rows)
- warm-start cache: dict `(gap, Q_total, eps_org) → DOF vector`
"""
mutable struct FitContext
    rows::Vector{BloomquistRow}
    sel::Vector{Int}                   # indices into rows[] selected for fit
    c_eq::Any                          # NamedTuple from solve_phosphate_equilibrium
    mesh_by_delta::Dict{Float64,Any}   # δ_round (μm) → mesh
    warm_by_key::Dict{Tuple{Float64,Float64,Float64},Vector{Float64}}
    N_mesh::Int
    stretch::Float64
    V_cathode_warm::Dict{Tuple{Float64,Float64,Float64},Float64}
end

"""
    build_context(rows, sel; N_mesh=100, stretch=10.0, warm_init=nothing)

Compute derived quantities (gap_m, j_target_A_m2, delta_lev_m, R_series_Ohm_m2)
for every row, build the bulk equilibrium once, and prepare a per-δ mesh
cache. `warm_init` may be a previously-computed Stage 3 cache mapping
(gap, Q_total, ε_org) → DOF vector; otherwise the context starts empty and
each (gap, Q, ε) tuple gets warm-started by Stage 4 on first encounter.
"""
function build_context(rows_raw::Vector{BloomquistRow}, sel::Vector{Int};
                       N_mesh::Int = 100,
                       stretch::Float64 = 10.0,
                       warm_init = nothing)
    c_eq = Chemistry.solve_phosphate_equilibrium()

    rows = similar(rows_raw)
    for (i, r) in pairs(rows_raw)
        gap_m   = r.gap_mm * 1.0e-3
        Qt_m3s  = ml_min_to_m3_s(r.Q_total_mL_min)
        δ_lev   = delta_leveque(gap_m, Qt_m3s)
        Rs      = R_series(gap_m, δ_lev, c_eq, r.phi_AN;
                           R_contact = R_CONTACT_DEFAULT)
        rows[i] = BloomquistRow(r.table, r.gap_mm, r.Q_total_mL_min,
                                r.j_mA_cm2, r.phi_AN,
                                r.Q_aq_mL_min, r.Q_org_mL_min,
                                r.We_aq, r.We_org,
                                r.FE_ADN_pct, r.FE_TCH_pct, r.FE_PN_pct,
                                r.PR_ADN_kg_cm2_h, r.EP_ADN_kg_kWh,
                                gap_m, r.j_mA_cm2 * 10.0,
                                δ_lev, Rs)
    end

    mesh_by_delta = Dict{Float64,Any}()
    for r in rows[sel]
        key = round(r.delta_lev_m * 1e6; digits = 1)   # μm rounded
        if !haskey(mesh_by_delta, key)
            mesh_by_delta[key] = make_mesh(N_mesh, r.delta_lev_m;
                                           stretch = stretch)
        end
    end

    warm_by_key = warm_init === nothing ?
        Dict{Tuple{Float64,Float64,Float64},Vector{Float64}}() :
        warm_init
    V_warm = Dict{Tuple{Float64,Float64,Float64},Float64}()

    return FitContext(rows, sel, c_eq, mesh_by_delta, warm_by_key,
                      N_mesh, stretch, V_warm)
end

@inline function _row_key(r::BloomquistRow)
    return (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
end

@inline function _mesh_key(r::BloomquistRow)
    return round(r.delta_lev_m * 1e6; digits = 1)
end

# ---------- Residual builder ----------
"""
    residuals!(F, theta, ctx; verbose=false) -> nrow_failed::Int

Fill `F` (length 2 · length(ctx.sel)) in-place with stacked
`(FE_ADN_model − FE_ADN_obs, FE_PN_model − FE_PN_obs)` residuals in
percentage points. Rows where the fixed-j solve fails to converge get a
large penalty (1e3 pp) on both residual entries and contribute to the
returned `nrow_failed` count. Updates the warm-start cache in-place.
"""
function residuals!(F::AbstractVector{Float64},
                    theta::AbstractVector{Float64},
                    ctx::FitContext;
                    verbose::Bool = false)
    @assert length(F) == 2 * length(ctx.sel)
    j0, ac = theta_to_physical(theta)

    nfail = 0
    @inbounds for (n, idx) in pairs(ctx.sel)
        r       = ctx.rows[idx]
        mesh_k  = _mesh_key(r)
        mesh    = ctx.mesh_by_delta[mesh_k]
        wkey    = _row_key(r)
        u_warm  = if haskey(ctx.warm_by_key, wkey)
                      ctx.warm_by_key[wkey]
                  else
                      # Cold start: zero-current equilibrium IC at this ε_org
                      make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN)
                  end

        result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                            mesh, u_warm, ctx.c_eq;
                            j0 = j0, alpha_c = ac,
                            verbose = false)

        if result.converged
            ctx.warm_by_key[wkey]      = result.state
            ctx.V_cathode_warm[wkey]   = result.V_cathode
            F[2n - 1] = result.FE_ADN_pct - r.FE_ADN_pct
            F[2n    ] = result.FE_PN_pct  - r.FE_PN_pct
        else
            nfail += 1
            verbose && @warn "fixed-j FAIL" row=idx note=result.note j_target=r.j_target_A_m2 ε=r.phi_AN gap=r.gap_mm
            F[2n - 1] = 1.0e3
            F[2n    ] = 1.0e3
        end
    end
    return nfail
end

"Sum-of-squared-residuals (pp²)."
function loss(theta::AbstractVector{Float64}, ctx::FitContext)
    F = zeros(2 * length(ctx.sel))
    nfail = residuals!(F, theta, ctx)
    return (loss = dot(F, F), nfail = nfail)
end

# ---------- Minimal pure-Julia Levenberg-Marquardt ----------
struct LMResult
    converged::Bool
    theta::Vector{Float64}
    loss::Float64
    iter::Int
    nfail_rows::Int
    history::Vector{NamedTuple{(:iter, :loss, :lambda, :nfail), Tuple{Int,Float64,Float64,Int}}}
    note::String
end

"""
    lm_fit(theta0, ctx; lb=THETA_LB, ub=THETA_UB,
           lambda0=1e-2, max_iter=80, tol_rel=1e-4,
           fd_step=1e-3, verbose=true) -> LMResult

Minimal damped Gauss–Newton (Levenberg–Marquardt) on the kinetics-fit problem.
Finite-difference Jacobian (forward differences, step `fd_step` in θ-units)
because each residual evaluation is itself a Newton solve — autodiff through
the inner solver isn't tractable here. With N_THETA = 6, J-build costs 6
loss evaluations per LM step.

Box bounds enforced by simple projection at the end of each accepted step.
"""
function lm_fit(theta0::AbstractVector{Float64},
                ctx::FitContext;
                lb::Vector{Float64}      = THETA_LB,
                ub::Vector{Float64}      = THETA_UB,
                lambda0::Float64         = 1.0e-2,
                lambda_up::Float64       = 4.0,
                lambda_down::Float64     = 0.5,
                max_iter::Int            = 80,
                tol_rel::Float64         = 1.0e-4,
                fd_step::Float64         = 1.0e-3,
                verbose::Bool            = true)

    @assert length(theta0) == N_THETA == length(lb) == length(ub)
    theta = clamp.(collect(theta0), lb, ub)
    nres  = 2 * length(ctx.sel)

    F = zeros(nres)
    nfail = residuals!(F, theta, ctx; verbose = false)
    loss_now = dot(F, F)
    history = NamedTuple{(:iter, :loss, :lambda, :nfail),
                         Tuple{Int,Float64,Float64,Int}}[
        (iter = 0, loss = loss_now, lambda = lambda0, nfail = nfail)
    ]

    if verbose
        @printf("[LM] iter %3d   loss=%.4e   λ=%.2e   nfail=%d\n",
                0, loss_now, lambda0, nfail)
    end

    lambda = lambda0
    J = zeros(nres, N_THETA)
    F_pert = zeros(nres)

    for it in 1:max_iter
        # ---- Forward-difference Jacobian ----
        for k in 1:N_THETA
            dθ = max(fd_step, fd_step * abs(theta[k]))
            θk_save = theta[k]
            theta[k] = clamp(θk_save + dθ, lb[k], ub[k])
            residuals!(F_pert, theta, ctx; verbose = false)
            J[:, k] = (F_pert .- F) ./ (theta[k] - θk_save + eps())
            theta[k] = θk_save
        end

        # ---- Damped Gauss–Newton step ----
        JtJ = J' * J
        Jtr = J' * F
        # Marquardt scaling: λ on the diagonal of JtJ
        H   = JtJ + lambda * Diagonal(diag(JtJ) .+ 1e-12)
        local Δ
        try
            Δ = -(H \ Jtr)
        catch err
            return LMResult(false, theta, loss_now, it, nfail, history,
                            "Linear solve failed: $(err)")
        end

        theta_trial = clamp.(theta .+ Δ, lb, ub)
        F_trial = zeros(nres)
        nfail_trial = residuals!(F_trial, theta_trial, ctx; verbose = false)
        loss_trial = dot(F_trial, F_trial)

        if loss_trial < loss_now
            # accept
            rel_drop = (loss_now - loss_trial) / max(loss_now, 1e-12)
            theta   = theta_trial
            F      .= F_trial
            loss_now = loss_trial
            nfail    = nfail_trial
            lambda   = max(lambda * lambda_down, 1e-8)
            push!(history, (iter = it, loss = loss_now, lambda = lambda,
                            nfail = nfail))
            if verbose
                @printf("[LM] iter %3d   loss=%.4e   λ=%.2e   nfail=%d   accept (Δ_rel=%.2e)\n",
                        it, loss_now, lambda, nfail, rel_drop)
            end
            if rel_drop < tol_rel
                return LMResult(true, theta, loss_now, it, nfail, history,
                                "Converged: relative loss drop < tol_rel")
            end
        else
            lambda = min(lambda * lambda_up, 1e8)
            push!(history, (iter = it, loss = loss_now, lambda = lambda,
                            nfail = nfail))
            if verbose
                @printf("[LM] iter %3d   loss=%.4e   λ=%.2e   nfail=%d   reject\n",
                        it, loss_now, lambda, nfail)
            end
            if lambda > 1e7
                return LMResult(false, theta, loss_now, it, nfail, history,
                                "Damping λ exceeded 1e7 — stuck")
            end
        end
    end

    return LMResult(false, theta, loss_now, max_iter, nfail, history,
                    "Reached max_iter without convergence")
end

end # module
