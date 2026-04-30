module FitKineticsSeq

# v7 §22 — Lees S1.11 sequential per-branch kinetics fit.
#
# Fits the four cathodic branches one at a time:
#     HER → ADN → PN → TCH
# At each sub-fit, the partial currents of all already-fit branches AND all
# not-yet-fit branches are clamped to their experimental values via
# Kinetics.FIXED_J.  Only the active branch's kinetics are varied; the loss is
# the sum of squared log10(j_branch + ε) residuals over Core rows.
#
# Why HER first: it is the only branch with no AN dependence (single-reactant
# Tafel), the analog of Lees' H₂ choice in the original CO2R workflow.
#
# Outer iteration: after one full HER→ADN→PN→TCH pass, optionally re-fit each
# branch with the others held at their newly-fitted values; repeat until
# parameter changes < tol_outer.
#
# Order of params per branch:
#   HER : (log10 j0_HER, alpha_c_HER)               2 dof
#   ADN : (log10 j0_ADN, alpha_c_ADN, n_ADN)        3 dof
#   PN  : (log10 j0_PN,  alpha_c_PN,  n_PN)         3 dof
#   TCH : (log10 j0_TCH, alpha_c_TCH, n_TCH)        3 dof

using Printf
using LinearAlgebra
using ..Params
using ..Mesh
using ..Chemistry
using ..Hydrodynamics
using ..CellVoltage
using ..Kinetics
using ..GalvContinuation
using ..FitKinetics: BloomquistRow, FitContext, build_context,
                     select_core, select_extended, select_holdout,
                     _row_key, _mesh_key
import ..FitKinetics

export SeqFitState, SeqFitState_from_v3,
       run_sequential_fit, branch_residuals!,
       J0_INIT_HER, AC_INIT_HER

# ------------------------------------------------------------
# Sub-fit parameter bounds and initial guesses
# ------------------------------------------------------------
# HER seed: v6.x converged values that the v3 fit was using as frozen.
const J0_INIT_HER = 2.666e-5      # A m⁻²
const AC_INIT_HER = 0.390

const HER_LB = [-7.0, 0.20]
const HER_UB = [-3.0, 0.70]
const HER_INIT = [log10(J0_INIT_HER), AC_INIT_HER]

# ADN/PN/TCH bounds — same as v3 single-branch slices (ADN n bottom = 1.0
# per Mathison JACS 2025; PN/TCH allowed down to 0.5 / 1.0 per v6.x history).
const ADN_LB = [-6.0, 0.30, 1.0];  ADN_UB = [-1.0, 0.70, 3.0]
const PN_LB  = [-6.0, 0.30, 0.5];  PN_UB  = [-1.0, 0.70, 2.0]
const TCH_LB = [-6.0, 0.30, 1.0];  TCH_UB = [-1.0, 0.70, 3.0]

# log10(j+EPS_J) regression: ε in A/m² (= 0.1 mA/cm²) avoids -Inf when j→0
const EPS_J_REG = 1.0   # A/m²

# ------------------------------------------------------------
# Fit-state container
# ------------------------------------------------------------
mutable struct SeqFitState
    # All four branches in physical units. NaN entries mean "not yet fit"
    # (the seed value is used wherever needed before its first sub-fit).
    j0_ADN::Float64;  ac_ADN::Float64;  n_ADN::Float64
    j0_PN ::Float64;  ac_PN ::Float64;  n_PN ::Float64
    j0_HER::Float64;  ac_HER::Float64
    j0_TCH::Float64;  ac_TCH::Float64;  n_TCH::Float64
    # Per-row warm starts: (gap, Q_total, eps_org) → (u, V, j_total).  Only
    # written on success. Cleared at the start of each outer iteration to keep
    # the cache consistent with the current parameter set.
    warm::Dict{Tuple{Float64,Float64,Float64},Tuple{Vector{Float64},Float64,Float64}}
end

function SeqFitState()
    SeqFitState(NaN, NaN, NaN,    # ADN
                NaN, NaN, NaN,    # PN
                J0_INIT_HER, AC_INIT_HER,  # HER (seeded)
                NaN, NaN, NaN,    # TCH
                Dict{Tuple{Float64,Float64,Float64},
                     Tuple{Vector{Float64},Float64,Float64}}())
end

"""
    SeqFitState_from_v3() -> SeqFitState

Seed all 11 kinetic parameters with the converged v3 fit (Stage 4v3 fitted_theta.txt).
Used as the warm-start initial guess for the sequential fit so that outer-1 only
needs to refine, not search. v3 results: HER frozen at v6.x, ADN/PN/TCH joint LM.
"""
function SeqFitState_from_v3()
    st = SeqFitState()
    # v3 converged values (from output/stage4v3/data/stage4a_fitted_theta.txt)
    st.j0_ADN = 8.454694e-3;  st.ac_ADN = 0.566594;  st.n_ADN = 1.000000
    st.j0_PN  = 1.616525e-3;  st.ac_PN  = 0.560814;  st.n_PN  = 1.636977
    st.j0_HER = 2.666e-5;     st.ac_HER = 0.390000     # frozen in v3, fittable here
    st.j0_TCH = 2.426481e-3;  st.ac_TCH = 0.563606;  st.n_TCH = 1.549037
    return st
end

# ------------------------------------------------------------
# Branch ordering: which slot is each branch in (j0/ac/n) tuples for
# tafel_currents / FIXED_J ?  (1=ADN, 2=PN, 3=HER, 4=TCH)
# ------------------------------------------------------------
const BRANCH_IDX = (ADN = 1, PN = 2, HER = 3, TCH = 4)

# ------------------------------------------------------------
# Build the (j0, ac, n) overrides for a given SeqFitState, with one branch
# replaced by the trial parameters.  Branches not yet fit fall back to seed
# defaults so the BV computation remains finite (their currents are clamped
# via FIXED_J anyway, so the BV value never enters the loss).
# ------------------------------------------------------------
function _build_kin_override(st::SeqFitState, branch::Symbol,
                              theta_branch::AbstractVector{Float64})
    # Defaults (also seed values for unfit branches) — chosen to not blow up BV.
    j0_def  = (10^-3, 10^-3, J0_INIT_HER, 10^-3)
    ac_def  = (0.5,   0.5,   AC_INIT_HER, 0.5)
    n_def   = (2.0,   1.0,   3.0)        # n is 3-tuple (ADN, PN, TCH)

    # Override with current state where known
    j0 = collect(j0_def);  ac = collect(ac_def);  n = collect(n_def)

    isnan(st.j0_ADN) || (j0[1] = st.j0_ADN; ac[1] = st.ac_ADN; n[1] = st.n_ADN)
    isnan(st.j0_PN ) || (j0[2] = st.j0_PN ; ac[2] = st.ac_PN ; n[2] = st.n_PN )
    j0[3] = st.j0_HER;  ac[3] = st.ac_HER
    isnan(st.j0_TCH) || (j0[4] = st.j0_TCH; ac[4] = st.ac_TCH; n[3] = st.n_TCH)

    # Active branch overrides
    if branch === :HER
        j0[3] = 10.0 ^ theta_branch[1]
        ac[3] =        theta_branch[2]
    elseif branch === :ADN
        j0[1] = 10.0 ^ theta_branch[1]
        ac[1] =        theta_branch[2]
        n[1]  =        theta_branch[3]
    elseif branch === :PN
        j0[2] = 10.0 ^ theta_branch[1]
        ac[2] =        theta_branch[2]
        n[2]  =        theta_branch[3]
    elseif branch === :TCH
        j0[4] = 10.0 ^ theta_branch[1]
        ac[4] =        theta_branch[2]
        n[3]  =        theta_branch[3]
    else
        error("unknown branch: $branch")
    end

    return ((j0[1], j0[2], j0[3], j0[4]),
            (ac[1], ac[2], ac[3], ac[4]),
            (n[1], n[2], n[3]))
end

# ------------------------------------------------------------
# Per-row partial currents from observation (Lees-clamped values).
# FE_HER is back-derived from closure: 100 - (FE_ADN + FE_PN + FE_TCH).
# Returns NaN for the active branch (we don't clamp it!).
# ------------------------------------------------------------
function _observed_partial_currents(r::BloomquistRow, active::Symbol)
    j_total = r.j_target_A_m2
    fe_her  = max(0.0, 100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct))
    j_ADN   = j_total * r.FE_ADN_pct / 100.0
    j_PN    = j_total * r.FE_PN_pct  / 100.0
    j_HER   = j_total * fe_her       / 100.0
    j_TCH   = j_total * r.FE_TCH_pct / 100.0
    # Clamp all branches except the active one
    out = (active === :ADN ? NaN : j_ADN,
           active === :PN  ? NaN : j_PN ,
           active === :HER ? NaN : j_HER,
           active === :TCH ? NaN : j_TCH)
    return out, j_total
end

# ------------------------------------------------------------
# Branch residuals: log10(j_branch_model + ε) − log10(j_branch_obs + ε)
# F has length = length(ctx.sel).  Returns nfail.
# ------------------------------------------------------------
"""
    branch_residuals!(F, theta, st, ctx, branch; verbose=false)
        -> (nfail::Int, fail_mask::BitVector)

Fill `F` with log10(j+ε) residuals for the active branch.

**Option 1 (current implementation): un-clamped sequential fit.**
The active branch's parameters come from `theta_branch`; all other branches
use whatever is currently committed in `st` (NaN sentinels fall back to
seeds). FIXED_J clamping is NOT used — the model's BV computes all four
partial currents from the kinetic params, j_total = j_target is enforced by
the bordered galvanostatic Newton, and the active branch's prediction is a
real function of theta_branch (so the gradient is non-zero).

This is "v3 joint LM, run one branch at a time" rather than true Lees S1.11
clamping, because clamping all-but-one branches AND enforcing j_total
makes the active-branch residual zero by closure (no fit information).

Failed rows: `F[i] = 0.0`, `fail_mask[i] = true` — LM drops them from J.
"""
function branch_residuals!(F::AbstractVector{Float64},
                            theta_branch::AbstractVector{Float64},
                            st::SeqFitState,
                            ctx::FitContext,
                            branch::Symbol;
                            verbose::Bool = false)
    @assert length(F) == length(ctx.sel)
    nfail = 0
    fail_mask = falses(length(ctx.sel))
    j0_use, ac_use, n_use = _build_kin_override(st, branch, theta_branch)

    @inbounds for (i, idx) in pairs(ctx.sel)
        r       = ctx.rows[idx]
        mesh_k  = _mesh_key(r)
        mesh    = ctx.mesh_by_delta[mesh_k]
        wkey    = _row_key(r)
        j_total = r.j_target_A_m2

        # Pull warm hint if we have one for this row
        has_warm = haskey(st.warm, wkey)
        warm_u = has_warm ? st.warm[wkey][1] : nothing
        warm_V = has_warm ? st.warm[wkey][2] : NaN
        warm_j = has_warm ? st.warm[wkey][3] : NaN

        # Un-clamped solve (Option 1) — no fixed_j passed
        result = solve_at_j_cont(j_total, r.phi_AN, mesh, ctx.c_eq;
                                 j0       = j0_use,
                                 alpha_c  = ac_use,
                                 n_orders = n_use,
                                 warm_u   = warm_u,
                                 warm_V   = warm_V,
                                 warm_j   = warm_j,
                                 verbose  = false)

        if result.converged
            st.warm[wkey] = (copy(result.state), result.V_cathode, result.j_total)
            j_branch_model = branch === :ADN ? result.j1 :
                             branch === :PN  ? result.j2 :
                             branch === :HER ? result.j3 :
                             branch === :TCH ? result.j4 :
                             error("bad branch")
            j_branch_obs = branch === :ADN ? j_total * r.FE_ADN_pct / 100.0 :
                           branch === :PN  ? j_total * r.FE_PN_pct  / 100.0 :
                           branch === :HER ? j_total * max(0.0,
                                                100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct)) / 100.0 :
                           branch === :TCH ? j_total * r.FE_TCH_pct / 100.0 :
                           error("bad branch")
            F[i] = log10(max(j_branch_model, 0.0) + EPS_J_REG) -
                   log10(max(j_branch_obs,   0.0) + EPS_J_REG)
        else
            nfail += 1
            fail_mask[i] = true
            verbose && @warn "seq solve FAIL" branch row=idx note=result.note
            F[i] = 0.0   # drop from loss; LM also drops this row from J
        end
    end
    return nfail, fail_mask
end

# ------------------------------------------------------------
# Per-branch LM (smaller than the v3 monolithic LM since dim is 2 or 3)
# ------------------------------------------------------------
struct BranchLMResult
    converged::Bool
    theta::Vector{Float64}
    loss::Float64
    iter::Int
    note::String
end

function _lm_branch(theta0::AbstractVector{Float64},
                    st::SeqFitState,
                    ctx::FitContext,
                    branch::Symbol;
                    lb::Vector{Float64},
                    ub::Vector{Float64},
                    lambda0::Float64 = 1.0e-2,
                    max_iter::Int    = 30,
                    tol_rel::Float64 = 5.0e-3,
                    fd_step::Float64 = 5.0e-3,
                    lambda_stuck::Float64 = 1.0e5,
                    verbose::Bool    = true)

    n_dof = length(theta0)
    @assert length(lb) == length(ub) == n_dof
    theta = clamp.(collect(theta0), lb, ub)
    nres  = length(ctx.sel)

    F = zeros(nres)
    nfail, fail_mask = branch_residuals!(F, theta, st, ctx, branch)
    loss_now = dot(F, F)

    if verbose
        @printf("[LM-%s] iter %3d   loss=%.4e   λ=%.2e   nfail=%d/%d\n",
                branch, 0, loss_now, lambda0, nfail, nres)
        flush(stdout)
    end

    lambda = lambda0
    J = zeros(nres, n_dof)
    F_pert = zeros(nres)

    for it in 1:max_iter
        # Forward-difference Jacobian.  For consistency with F we want the
        # SAME row mask (a row that solved at theta but fails at theta+dθ
        # would otherwise produce a huge spurious gradient).  We zero-out any
        # row that is masked at *either* the base θ or the perturbed θ.
        for k in 1:n_dof
            dθ = max(fd_step, fd_step * abs(theta[k]))
            θk_save = theta[k]
            theta[k] = clamp(θk_save + dθ, lb[k], ub[k])
            _, pert_mask = branch_residuals!(F_pert, theta, st, ctx, branch)
            J[:, k] = (F_pert .- F) ./ (theta[k] - θk_save + eps())
            # Drop rows that became unsolved at the perturbed θ
            for i in 1:nres
                if fail_mask[i] || pert_mask[i]
                    J[i, k] = 0.0
                end
            end
            theta[k] = θk_save
        end

        # Active-set projection
        Jtr = J' * F
        bnd_tol = 1e-8
        active = falses(n_dof)
        for k in 1:n_dof
            at_lb = theta[k] ≤ lb[k] + bnd_tol
            at_ub = theta[k] ≥ ub[k] - bnd_tol
            active[k] = (at_lb && Jtr[k] > 0.0) || (at_ub && Jtr[k] < 0.0)
        end
        free = findall(.!active)
        if isempty(free)
            return BranchLMResult(false, theta, loss_now, it,
                                  "All params at active bounds")
        end
        if verbose && any(active)
            @printf("[LM-%s] iter %3d   active bounds: %s\n", branch, it,
                    join(string.(findall(active)), ","))
        end

        # Damped GN on free params
        J_r   = J[:, free]
        JtJ_r = J_r' * J_r
        Jtr_r = Jtr[free]
        H_r   = JtJ_r + lambda * Diagonal(diag(JtJ_r) .+ 1e-12)
        Δ = zeros(n_dof)
        try
            Δ[free] = -(H_r \ Jtr_r)
        catch err
            return BranchLMResult(false, theta, loss_now, it,
                                  "linear solve failed: $err")
        end
        theta_trial = clamp.(theta .+ Δ, lb, ub)
        F_trial = zeros(nres)
        nfail_trial, mask_trial = branch_residuals!(F_trial, theta_trial, st, ctx, branch)
        loss_trial = dot(F_trial, F_trial)

        if loss_trial < loss_now
            rel_drop = (loss_now - loss_trial) / max(loss_now, 1e-12)
            theta = theta_trial
            F .= F_trial
            loss_now = loss_trial
            nfail = nfail_trial
            fail_mask .= mask_trial      # carry mask to next iter's J build
            lambda = max(lambda * 0.5, 1e-8)
            if verbose
                @printf("[LM-%s] iter %3d   loss=%.4e   λ=%.2e   nfail=%d   accept (Δ_rel=%.2e)\n",
                        branch, it, loss_now, lambda, nfail, rel_drop)
            end
            if rel_drop < tol_rel
                return BranchLMResult(true, theta, loss_now, it,
                                      "converged: rel drop < tol")
            end
        else
            lambda = min(lambda * 4.0, 1e6)
            if verbose
                @printf("[LM-%s] iter %3d   loss=%.4e   λ=%.2e   nfail=%d   reject\n",
                        branch, it, loss_now, lambda, nfail)
            end
            if lambda > lambda_stuck
                return BranchLMResult(false, theta, loss_now, it,
                                      "λ > $(lambda_stuck) — stuck")
            end
        end
    end
    return BranchLMResult(false, theta, loss_now, max_iter,
                          "max_iter reached")
end

# ------------------------------------------------------------
# Commit a converged branch result back to the state
# ------------------------------------------------------------
function _commit!(st::SeqFitState, branch::Symbol, theta::Vector{Float64})
    if branch === :HER
        st.j0_HER = 10^theta[1];  st.ac_HER = theta[2]
    elseif branch === :ADN
        st.j0_ADN = 10^theta[1];  st.ac_ADN = theta[2];  st.n_ADN = theta[3]
    elseif branch === :PN
        st.j0_PN  = 10^theta[1];  st.ac_PN  = theta[2];  st.n_PN  = theta[3]
    elseif branch === :TCH
        st.j0_TCH = 10^theta[1];  st.ac_TCH = theta[2];  st.n_TCH = theta[3]
    else
        error("unknown branch: $branch")
    end
end

# ------------------------------------------------------------
# Top-level driver: HER → ADN → PN → TCH, optional outer iterations.
# ------------------------------------------------------------
"""
    run_sequential_fit(ctx; max_outer=5, branch_max_iter=30, verbose=true,
                        warm_start::Union{Nothing,SeqFitState}=nothing,
                        lambda_stuck=1e5)
        -> (state, history)

Sequential per-branch fit (HER → ADN → PN → TCH) on the rows in `ctx.sel`.
Returns the SeqFitState **with the lowest total loss across outer iterations**
(not necessarily the last one — outer iteration can drift) and a Vector of
records.

`warm_start`: optional SeqFitState to seed all 11 params (e.g., v3 results).
`lambda_stuck`: per-branch LM termination threshold (passed to `_lm_branch`).
"""
function run_sequential_fit(ctx::FitContext;
                            max_outer::Int       = 5,
                            branch_max_iter::Int = 30,
                            tol_outer::Float64   = 5.0e-3,
                            warm_start::Union{Nothing,SeqFitState} = nothing,
                            lambda_stuck::Float64 = 1.0e5,
                            verbose::Bool        = true)
    st = warm_start === nothing ? SeqFitState() :
            SeqFitState(warm_start.j0_ADN, warm_start.ac_ADN, warm_start.n_ADN,
                        warm_start.j0_PN,  warm_start.ac_PN,  warm_start.n_PN,
                        warm_start.j0_HER, warm_start.ac_HER,
                        warm_start.j0_TCH, warm_start.ac_TCH, warm_start.n_TCH,
                        Dict{Tuple{Float64,Float64,Float64},
                             Tuple{Vector{Float64},Float64,Float64}}())
    history = NamedTuple{(:outer, :branch, :loss, :iter, :note),
                         Tuple{Int,Symbol,Float64,Int,String}}[]

    # ---- Best-loss-so-far tracking (improvement #3) ----
    # Snapshot the parameter set after every outer iteration; keep the one
    # whose summed branch losses are smallest. This protects against
    # outer-iteration drift to a worse fit.
    best_st         = deepcopy(st)
    best_total_loss = Inf
    best_outer_idx  = 0

    # Param snapshot for outer-iter convergence check
    snapshot_prev = (st.j0_HER, st.ac_HER,
                     st.j0_ADN, st.ac_ADN, st.n_ADN,
                     st.j0_PN,  st.ac_PN,  st.n_PN,
                     st.j0_TCH, st.ac_TCH, st.n_TCH)

    for outer in 1:max_outer
        verbose && @info "===== Outer iteration $outer / $max_outer ====="
        # Wipe warm cache between outer passes (params change ⇒ stale)
        empty!(st.warm)

        # ---- HER ----
        seed = isnan(st.j0_HER) ? HER_INIT : [log10(st.j0_HER), st.ac_HER]
        verbose && @info "Sub-fit: HER" seed
        rh = _lm_branch(seed, st, ctx, :HER;
                        lb = HER_LB, ub = HER_UB,
                        max_iter = branch_max_iter,
                        lambda_stuck = lambda_stuck,
                        verbose = verbose)
        _commit!(st, :HER, rh.theta)
        push!(history, (outer = outer, branch = :HER, loss = rh.loss,
                        iter = rh.iter, note = rh.note))

        # ---- ADN ----
        adn_init = isnan(st.j0_ADN) ? [-3.0, 0.55, 1.0] :
                                      [log10(st.j0_ADN), st.ac_ADN, st.n_ADN]
        verbose && @info "Sub-fit: ADN" seed=adn_init
        ra = _lm_branch(adn_init, st, ctx, :ADN;
                        lb = ADN_LB, ub = ADN_UB,
                        max_iter = branch_max_iter,
                        lambda_stuck = lambda_stuck,
                        verbose = verbose)
        _commit!(st, :ADN, ra.theta)
        push!(history, (outer = outer, branch = :ADN, loss = ra.loss,
                        iter = ra.iter, note = ra.note))

        # ---- PN ----
        pn_init = isnan(st.j0_PN) ? [-3.0, 0.55, 1.0] :
                                    [log10(st.j0_PN), st.ac_PN, st.n_PN]
        verbose && @info "Sub-fit: PN" seed=pn_init
        rp = _lm_branch(pn_init, st, ctx, :PN;
                        lb = PN_LB, ub = PN_UB,
                        max_iter = branch_max_iter,
                        lambda_stuck = lambda_stuck,
                        verbose = verbose)
        _commit!(st, :PN, rp.theta)
        push!(history, (outer = outer, branch = :PN, loss = rp.loss,
                        iter = rp.iter, note = rp.note))

        # ---- TCH ----
        tch_init = isnan(st.j0_TCH) ? [-3.0, 0.55, 1.5] :
                                      [log10(st.j0_TCH), st.ac_TCH, st.n_TCH]
        verbose && @info "Sub-fit: TCH" seed=tch_init
        rt = _lm_branch(tch_init, st, ctx, :TCH;
                        lb = TCH_LB, ub = TCH_UB,
                        max_iter = branch_max_iter,
                        lambda_stuck = lambda_stuck,
                        verbose = verbose)
        _commit!(st, :TCH, rt.theta)
        push!(history, (outer = outer, branch = :TCH, loss = rt.loss,
                        iter = rt.iter, note = rt.note))

        # Best-loss-so-far check (improvement #3)
        outer_total_loss = rh.loss + ra.loss + rp.loss + rt.loss
        verbose && @info "Outer total loss" outer outer_total_loss best_total_loss
        if outer_total_loss < best_total_loss
            best_total_loss = outer_total_loss
            best_st         = deepcopy(st)
            best_outer_idx  = outer
        end

        # Outer convergence check
        snapshot_now = (st.j0_HER, st.ac_HER,
                        st.j0_ADN, st.ac_ADN, st.n_ADN,
                        st.j0_PN,  st.ac_PN,  st.n_PN,
                        st.j0_TCH, st.ac_TCH, st.n_TCH)
        rel_change = maximum(abs.(collect(snapshot_now) .- collect(snapshot_prev)) ./
                             max.(abs.(collect(snapshot_prev)), 1e-8))
        verbose && @info "Outer change" rel_change
        if rel_change < tol_outer
            verbose && @info "Outer iteration converged (rel_change < $tol_outer)"
            break
        end
        snapshot_prev = snapshot_now
    end

    verbose && @info "Returning BEST-LOSS state" best_outer_idx best_total_loss
    return best_st, history
end

end # module
