module FitKinetics

# v7 §20 — kinetics-only fit on Bloomquist data with all transport frozen.
#
#   theta = (log10 j0_ADN, log10 j0_PN, log10 j0_TCH,
#            alpha_c_ADN, alpha_c_PN, alpha_c_TCH,
#            n_ADN, n_PN, n_TCH)
#
# HER (j0_3, alpha_c3) frozen at v6.x converged values and injected via
# J0_3_FROZEN / ALPHA_C3_FROZEN — not fit parameters.
# Residual vector: 3 · length(sel) with [FE_ADN, FE_PN, FE_TCH] per row.
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
       N_THETA, THETA_LB, THETA_UB, THETA0,
       J0_3_FROZEN, ALPHA_C3_FROZEN,
       SOLVER_BACKEND

# v7 §21 — switch between legacy V-walk solver (`:vwalk`) and the new
# galvanostatic continuation solver (`:cont`). The continuation solver is
# faster and avoids warm-start corruption under FD-Jacobian perturbations.
# Set via `FitKinetics.SOLVER_BACKEND[] = :cont` BEFORE building the context.
const SOLVER_BACKEND = Ref{Symbol}(:vwalk)

# ---------- Frozen HER parameters (v6.x converged, not fit in v7) ----------
const J0_3_FROZEN     = 2.666e-5   # A m⁻²
const ALPHA_C3_FROZEN = 0.390

# ---------- Parameter vector layout ----------
const N_THETA = 9
# theta is in *log10* space for j0,r and *linear* for α_c,r and n_r.
# Layout: (log10 j0_ADN, log10 j0_PN, log10 j0_TCH,
#          ac_ADN, ac_PN, ac_TCH, n_ADN, n_PN, n_TCH)
# HER slot is absent — frozen via J0_3_FROZEN / ALPHA_C3_FROZEN above.
const THETA_LB = [-6.0, -6.0, -6.0, 0.30, 0.30, 0.30, 1.0, 0.5, 1.0]
const THETA_UB = [-1.0, -1.0, -1.0, 0.70, 0.70, 0.70, 3.0, 2.0, 3.0]
const THETA0   = [-3.0, -3.0, -3.0, 0.50, 0.50, 0.50, 2.0, 1.0, 2.0]

"Convert internal theta vector → physical (j0_4tuple, ac_4tuple, n_3tuple)."
function theta_to_physical(theta::AbstractVector{Float64})
    # HER (slot 3) injected from frozen constants; not in theta.
    j0 = (10.0^theta[1], 10.0^theta[2], J0_3_FROZEN,     10.0^theta[3])
    ac = (theta[4],      theta[5],       ALPHA_C3_FROZEN, theta[6])
    n  = (theta[7], theta[8], theta[9])
    return j0, ac, n
end

"Inverse of theta_to_physical for diagnostics / initial guesses."
function physical_to_theta(j0::NTuple{4,Float64}, ac::NTuple{4,Float64},
                           n::NTuple{3,Float64} = (2.0, 1.0, 3.0))
    return [log10(j0[1]), log10(j0[2]), log10(j0[4]),
            ac[1], ac[2], ac[4], n[1], n[2], n[3]]
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

# ---------- On-disk cache helpers (v5-compatible filename) ----------
const CACHE_DIR_DEFAULT = joinpath(@__DIR__, "output", "cache")

"v5/v6 cache filename: s_eo<ε>_d<δ_μm>_V<V>.bin"
function _cache_filename(eps_org::Float64, delta_m::Float64, V::Float64)
    delta_um = round(Int, delta_m * 1e6)
    return Printf.@sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, delta_um, V)
end

"Load DOF vector from a v5/v6 cache file. Returns nothing if file missing."
function _load_cache_state(path::String)
    isfile(path) || return nothing
    open(path, "r") do f
        n = read(f, Int64)
        u = Vector{Float64}(undef, n)
        read!(f, u)
        return u
    end
end

"""
    build_context(rows, sel; N_mesh=100, stretch=10.0,
                  warm_init=nothing, cache_dir=CACHE_DIR_DEFAULT,
                  V_warm=-1.0)

Compute derived quantities (gap_m, j_target_A_m2, delta_lev_m, R_series_Ohm_m2)
for every row, build the bulk equilibrium once, and prepare a per-δ mesh
cache. The warm-start cache is populated in this priority order:

  1. `warm_init` dict (if passed) — overrides everything else.
  2. On-disk Stage 3 cache files at `cache_dir`/`s_eo<ε>_d<δ_μm>_V<V_warm>.bin`,
     loaded for each unique (gap, Q_total, ε_org) key the selection touches.
  3. Cold start via `make_initial_guess` inside `residuals!` for any key still
     uncached at first encounter.
"""
function build_context(rows_raw::Vector{BloomquistRow}, sel::Vector{Int};
                       N_mesh::Int = 100,
                       stretch::Float64 = 10.0,
                       warm_init = nothing,
                       cache_dir::String = CACHE_DIR_DEFAULT,
                       V_warm::Float64 = -1.0)
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

    # If no in-memory warm-start dict was passed, try the on-disk Stage 3 cache.
    if warm_init === nothing && isdir(cache_dir)
        n_loaded = 0
        for r in rows[sel]
            wkey = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
            haskey(warm_by_key, wkey) && continue
            fname = _cache_filename(r.phi_AN, r.delta_lev_m, V_warm)
            u = _load_cache_state(joinpath(cache_dir, fname))
            if u !== nothing && length(u) == 10 * N_mesh
                warm_by_key[wkey] = u
                n_loaded += 1
            end
        end
        if n_loaded > 0
            @info "FitContext: loaded $(n_loaded) warm-start states from $(cache_dir)"
        end
    end

    V_cathode_warm = Dict{Tuple{Float64,Float64,Float64},Float64}()

    return FitContext(rows, sel, c_eq, mesh_by_delta, warm_by_key,
                      N_mesh, stretch, V_cathode_warm)
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
    @assert length(F) == 3 * length(ctx.sel)
    j0, ac, n_orders = theta_to_physical(theta)
    backend = SOLVER_BACKEND[]

    nfail = 0
    @inbounds for (n, idx) in pairs(ctx.sel)
        r       = ctx.rows[idx]
        mesh_k  = _mesh_key(r)
        mesh    = ctx.mesh_by_delta[mesh_k]
        wkey    = _row_key(r)
        has_warm   = haskey(ctx.warm_by_key, wkey)
        has_warm_V = haskey(ctx.V_cathode_warm, wkey)
        u_warm  = has_warm ? ctx.warm_by_key[wkey] :
                             make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN)

        if backend === :cont
            # New galvanostatic continuation solver. Warm-start hint is
            # (u, V, j_total) from the previous successful solve at this row.
            warm_j = has_warm_V ? r.j_target_A_m2 : NaN   # j matches target
            # Resolve the continuation solver via Base.invokelatest to avoid
            # the FitKinetics-loaded-before-GalvContinuation load-order issue.
            result = Base.invokelatest(
                getfield(parentmodule(@__MODULE__).GalvContinuation,
                         :solve_at_j_cont),
                r.j_target_A_m2, r.phi_AN, mesh, ctx.c_eq;
                j0 = j0, alpha_c = ac, n_orders = n_orders,
                warm_u = has_warm_V ? u_warm : nothing,
                warm_V = has_warm_V ? ctx.V_cathode_warm[wkey] : NaN,
                warm_j = warm_j,
                verbose = false)
        else
            # Legacy V-walk + bisection solver (default).
            result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                                mesh, u_warm, ctx.c_eq;
                                j0 = j0, alpha_c = ac, n_orders = n_orders,
                                V_hi  = has_warm_V ? -0.8 : -0.3,
                                V_hint = has_warm_V ? ctx.V_cathode_warm[wkey] : NaN,
                                u_hint = has_warm_V ? u_warm : nothing,
                                verbose = false)
        end

        if result.converged
            ctx.warm_by_key[wkey]    = result.state
            ctx.V_cathode_warm[wkey] = result.V_cathode
            F[3n - 2] = result.FE_ADN_pct - r.FE_ADN_pct
            F[3n - 1] = result.FE_PN_pct  - r.FE_PN_pct
            F[3n    ] = result.FE_TCH_pct - r.FE_TCH_pct
        else
            nfail += 1
            verbose && @warn "fixed-j FAIL" row=idx note=result.note j_target=r.j_target_A_m2 ε=r.phi_AN gap=r.gap_mm
            F[3n - 2] = 1000.0
            F[3n - 1] = 1000.0
            F[3n    ] = 1000.0
        end
    end
    return nfail
end

"Sum-of-squared-residuals (pp²)."
function loss(theta::AbstractVector{Float64}, ctx::FitContext)
    F = zeros(3 * length(ctx.sel))
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
    F::Vector{Float64}   # residual vector at best theta (avoids re-evaluation after LM)
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
                lambda_stuck::Float64    = 1.0e3,    # was 1e7 — loosened for v6.x test-drives
                max_iter::Int            = 80,
                tol_rel::Float64         = 1.0e-2,   # was 1e-4 — loosened for v6.x test-drives
                fd_step::Float64         = 1.0e-3,
                # ---- L2 regularization toward a prior θ (added v7 §22) ----
                # When `theta_prior` is non-`nothing`, the residual vector is
                # augmented with `sqrt(prior_weight[k]) · (θ[k] - prior[k])`
                # for each parameter k. Loss = ||F_data||² + Σ w_k (θ_k - p_k)².
                # Useful for "warm-anchored" fits that want post-fix corrections
                # without abandoning a known-good prior.
                theta_prior::Union{Nothing,AbstractVector{Float64}} = nothing,
                prior_weight::Union{Float64,AbstractVector{Float64}} = 0.0,
                verbose::Bool            = true)

    @assert length(theta0) == N_THETA == length(lb) == length(ub)
    theta = clamp.(collect(theta0), lb, ub)

    # Build prior-weight vector (broadcast scalar; allow per-param weights).
    use_prior = theta_prior !== nothing && !(prior_weight isa Number && prior_weight == 0.0)
    pw_vec = use_prior ?
                (prior_weight isa Number ? fill(Float64(prior_weight), N_THETA) :
                                            collect(Float64.(prior_weight))) :
                zeros(Float64, N_THETA)
    pr_vec = use_prior ? collect(Float64.(theta_prior)) : zeros(Float64, N_THETA)
    @assert length(pw_vec) == N_THETA
    @assert length(pr_vec) == N_THETA

    n_data = 3 * length(ctx.sel)
    nres   = use_prior ? n_data + N_THETA : n_data

    # `residuals_full!` writes data residuals to F[1:n_data] and (if regularizing)
    # the regularization tail F[n_data+1 : n_data+N_THETA].
    function residuals_full!(F_full::AbstractVector{Float64},
                              th::AbstractVector{Float64})
        nf = residuals!(view(F_full, 1:n_data), th, ctx; verbose = false)
        if use_prior
            for k in 1:N_THETA
                F_full[n_data + k] = sqrt(pw_vec[k]) * (th[k] - pr_vec[k])
            end
        end
        return nf
    end

    F = zeros(nres)
    nfail = residuals_full!(F, theta)
    loss_now = dot(F, F)
    history = NamedTuple{(:iter, :loss, :lambda, :nfail),
                         Tuple{Int,Float64,Float64,Int}}[
        (iter = 0, loss = loss_now, lambda = lambda0, nfail = nfail)
    ]

    if verbose
        @printf("[LM] iter %3d   loss=%.4e   λ=%.2e   nfail=%d\n",
                0, loss_now, lambda0, nfail)
        flush(stdout)
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
            residuals_full!(F_pert, theta)
            J[:, k] = (F_pert .- F) ./ (theta[k] - θk_save + eps())
            theta[k] = θk_save
        end

        # ---- Active-set detection: freeze parameters at bounds where the ----
        # ---- gradient pushes further into the bound (KKT condition).      ----
        Jtr = J' * F
        bnd_tol = 1e-8
        active = falses(N_THETA)
        for k in 1:N_THETA
            at_lb = theta[k] ≤ lb[k] + bnd_tol
            at_ub = theta[k] ≥ ub[k] - bnd_tol
            # gradient component Jtr[k] = ∂loss/∂theta[k]; descent is -Jtr[k].
            # Frozen when descent direction would violate the bound.
            active[k] = (at_lb && Jtr[k] > 0.0) || (at_ub && Jtr[k] < 0.0)
        end
        free = findall(.!active)

        if isempty(free)
            return LMResult(false, theta, loss_now, it, nfail, history,
                            "All parameters at active bounds — KKT satisfied", copy(F))
        end

        if verbose && any(active)
            @printf("[LM] iter %3d   active bounds: %s\n", it,
                    join(string.(findall(active)), ","))
        end

        # ---- Damped Gauss–Newton step on free parameters only ----
        J_r   = J[:, free]
        JtJ_r = J_r' * J_r
        Jtr_r = Jtr[free]
        H_r   = JtJ_r + lambda * Diagonal(diag(JtJ_r) .+ 1e-12)
        Δ = zeros(N_THETA)
        try
            Δ[free] = -(H_r \ Jtr_r)
        catch err
            return LMResult(false, theta, loss_now, it, nfail, history,
                            "Linear solve failed: $(err)", copy(F))
        end

        theta_trial = clamp.(theta .+ Δ, lb, ub)
        F_trial = zeros(nres)
        nfail_trial = residuals_full!(F_trial, theta_trial)
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
                                "Converged: relative loss drop < tol_rel", copy(F))
            end
        else
            lambda = min(lambda * lambda_up, lambda_stuck * 10)
            push!(history, (iter = it, loss = loss_now, lambda = lambda,
                            nfail = nfail))
            if verbose
                @printf("[LM] iter %3d   loss=%.4e   λ=%.2e   nfail=%d   reject\n",
                        it, loss_now, lambda, nfail)
            end
            if lambda > lambda_stuck
                return LMResult(false, theta, loss_now, it, nfail, history,
                                "Damping λ exceeded lambda_stuck=$(lambda_stuck) — stuck (returning best so far)",
                                copy(F))
            end
        end
    end

    return LMResult(false, theta, loss_now, max_iter, nfail, history,
                    "Reached max_iter without convergence", copy(F))
end

end # module
