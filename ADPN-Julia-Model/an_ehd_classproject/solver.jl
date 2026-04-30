module Solver

using SparseArrays
using LinearAlgebra
using ForwardDiff

using ..Assembly

export build_banded_pattern,
       banded_fd_jacobian!,
       newton_solve!,
       newton_galvanostatic!,
       newton_continuation,
       newton_continuation_logj

# ------------------------------------------------------------
# Build sparsity pattern: block-tridiagonal, block size 10.
# Cell-major DOF layout ⇒ coupling of cell ix to cells ix−1, ix, ix+1
# ⇒ half-bandwidth b = 2·blocksize − 1 = 19.
# Column-grouped FD requires 2b+1 = 39 perturbation vectors.
# ------------------------------------------------------------
const JAC_BLOCK  = 10
const JAC_HALFBW = 2 * JAC_BLOCK - 1   # = 19

function build_banded_pattern(n::Int, b::Int = JAC_HALFBW)
    I = Int[]; J = Int[]; V = Float64[]
    for j in 1:n
        for i in max(1, j - b):min(n, j + b)
            push!(I, i); push!(J, j); push!(V, 0.0)
        end
    end
    return sparse(I, J, V, n, n)
end

# ------------------------------------------------------------
# Banded FD Jacobian with 19-color column grouping (§10.2).
# ------------------------------------------------------------
"""
    banded_fd_jacobian!(J, residual!, u, F0, b=9)

Overwrites J (SparseMatrixCSC with pre-built banded pattern) in place via
column-grouped forward FD. Requires 2b+1 = 19 residual evaluations.

Perturbations:
  log-concentration DOFs   ε = 1e-7 · max(|u_j|, 1)
  φ_l DOFs                 ε = 1e-5 V
"""
function banded_fd_jacobian!(J::SparseMatrixCSC{Float64,Int},
                             residual!::F, u::Vector{Float64},
                             F0::Vector{Float64}; b::Int = JAC_HALFBW) where {F<:Function}

    n        = length(u)
    n_colors = 2 * b + 1   # 19
    F_pert   = similar(F0)
    du       = zeros(n)

    # Zero out all non-zero entries (keep pattern)
    fill!(J.nzval, 0.0)

    for color in 1:n_colors
        # ---- build perturbation vector for this color ----
        fill!(du, 0.0)
        for j in color:n_colors:n
            is_phi = (j % 10 == 0)
            if is_phi
                du[j] = 1.0e-5
            else
                du[j] = 1.0e-7 * max(abs(u[j]), 1.0)
            end
        end

        # ---- evaluate residual at u + du ----
        u_pert = u .+ du
        residual!(F_pert, u_pert)

        # ---- scatter dF/du_j into J columns for this color ----
        for j in color:n_colors:n
            eps_j = du[j]
            eps_j == 0.0 && continue
            # rows affected: max(1, j-b) .. min(n, j+b)
            for i in max(1, j - b):min(n, j + b)
                J[i, j] = (F_pert[i] - F0[i]) / eps_j
            end
        end
    end
    return J
end

# ------------------------------------------------------------
# Newton solver with per-DOF step clamping (§10.3).
# ------------------------------------------------------------
"""
    newton_solve!(u, residual!; tol=1e-5, max_iter=40, jacobian_mode=:fd, verbose=false)

Newton with small diagonal Tikhonov `(J + λI) du = -F`, strict L2 descent +
backtracking, and optional FD or ForwardDiff AD Jacobian.

Fixed `λ = 1e-10` keeps the Newton direction physical in the near-singular
regime while preventing exact-singular linear solves. This is empirically
the branch that produces the physical Stage 1 solutions.

Step clamping: `Δu_log ≤ 5.0`, `Δφ_ℓ ≤ 0.015 V` per iteration.

`jacobian_mode`:
- `:fd` → banded column-grouped FD (35 perturbations, fast, cond-limited; default)
- `:ad` → ForwardDiff dense Jacobian (exact to machine precision, ~4× slower
  per Jacobian; use only for diagnostics or regions where FD fails)
"""
function newton_solve!(u::Vector{Float64}, residual!::F;
                       tol::Float64 = 1.0e-5,
                       max_iter::Int = 40,
                       max_backtrack::Int = 10,
                       lambda_fix::Float64 = 1.0e-10,
                       jacobian_mode::Symbol = :fd,
                       verbose::Bool = false) where {F<:Function}

    n       = length(u)
    F0      = zeros(n)
    Ftry    = zeros(n)
    N_cells = n ÷ 10   # v7: 10 DOF per cell (9 species + φ_l)

    Jsp     = build_banded_pattern(n, JAC_HALFBW)
    Jdense  = jacobian_mode === :ad ? zeros(n, n) : nothing
    ad_cfg  = jacobian_mode === :ad ? ForwardDiff.JacobianConfig(
                  residual!, zeros(n), u, ForwardDiff.Chunk{12}()) : nothing

    for iter in 1:max_iter
        residual!(F0, u)
        normF   = maximum(abs.(F0))
        normF_2 = sqrt(sum(abs2, F0))
        verbose && @info "Newton" iter normF normF_2
        if !isfinite(normF)
            return (converged = false, iter = iter, normF = normF)
        end
        if normF < tol
            return (converged = true, iter = iter, normF = normF)
        end

        # --- Build Jacobian ---
        if jacobian_mode === :ad
            ForwardDiff.jacobian!(Jdense, residual!, zeros(n), u, ad_cfg)
            # Add fixed Tikhonov damping to the diagonal.
            @inbounds for i in 1:n; Jdense[i, i] += lambda_fix; end
        else
            banded_fd_jacobian!(Jsp, residual!, u, F0)
            @inbounds for i in 1:n; Jsp[i, i] += lambda_fix; end
        end

        # --- Solve (J + λI) du = -F ---
        local du
        try
            if jacobian_mode === :ad
                du = -(Jdense \ F0)
            else
                du = -(Jsp \ F0)
            end
        catch e
            if e isa LinearAlgebra.SingularException ||
               (e isa LinearAlgebra.LAPACKException) ||
               occursin("Singular", sprint(showerror, e))
                verbose && @warn "Jacobian singular — aborting Newton" iter normF exception=e
                return (converged = false, iter = iter, normF = normF)
            else
                rethrow(e)
            end
        end
        if !all(isfinite, du)
            verbose && @warn "Newton step non-finite — aborting" iter normF
            return (converged = false, iter = iter, normF = normF)
        end

        # --- per-DOF-type step clamping (absolute bounds, §10.3) ---
        @inbounds for ix in 1:N_cells
            base = 10 * (ix - 1)
            for k in 1:9           # all 9 species (log-concentration DOFs)
                idx = base + k
                du[idx] = clamp(du[idx], -5.0, 5.0)
            end
            du[base + 10] = clamp(du[base + 10], -0.015, 0.015)  # φ_l
        end

        # --- Backtracking line search, strict L2 descent ---
        α = 1.0
        accepted = false
        for _ in 0:max_backtrack
            u_try = u .+ α .* du
            residual!(Ftry, u_try)
            normF_try_2 = sqrt(sum(abs2, Ftry))
            if isfinite(normF_try_2) && normF_try_2 < normF_2
                u .= u_try
                accepted = true
                break
            end
            α *= 0.5
        end
        if !accepted
            verbose && @warn "Newton line search failed — aborting" iter normF
            return (converged = false, iter = iter, normF = normF)
        end
    end

    residual!(F0, u)
    return (converged = false, iter = max_iter, normF = maximum(abs.(F0)))
end

# ------------------------------------------------------------
# Galvanostatic Newton — bordered (n+1)×(n+1) system (v7 §20.6)
# ------------------------------------------------------------
"""
    newton_galvanostatic!(u, V0, j_target, build_residual, j_total_fn; ...)

Newton solve for the galvanostatic BVP: find (u, V) such that
  F_phys(u, V) = 0   (the 10N-DOF physics residual)
  j_total(u, V) = j_target

The (n+1)×(n+1) bordered system is solved efficiently via the Schur complement:
only 2 banded LU solves are needed per Newton iteration instead of a dense O(n³) LU.

`build_residual(V)` → closure `(F, x) -> full_residual!(...)` at voltage V.
`j_total_fn(u, V)` → scalar j_total [A m⁻²].

DOF layout: 10 per cell (species 1–9 = log-concentrations, DOF 10 = φ_l).
Cell 1 is the electrode face (where Tafel BCs live).

Returns NamedTuple (converged, iter, normF, u, V).
"""
function newton_galvanostatic!(u::Vector{Float64},
                               V0::Float64,
                               j_target::Float64,
                               build_residual::BF,
                               j_total_fn::JF;
                               tol_phys::Float64   = 1.0e-5,
                               tol_rel_j::Float64  = 1.0e-3,
                               max_iter::Int       = 60,
                               max_backtrack::Int  = 12,
                               lambda_fix::Float64 = 1.0e-10,
                               verbose::Bool       = false) where {BF<:Function, JF<:Function}

    n       = length(u)
    N_cells = n ÷ 10          # v7: 10 DOF per cell
    V       = V0
    eps_V   = 1.0e-5

    Fphys = zeros(n)
    Fpert = zeros(n)
    Jsp   = build_banded_pattern(n, JAC_HALFBW)

    for iter in 1:max_iter
        residual! = build_residual(V)
        residual!(Fphys, u)
        j_now   = j_total_fn(u, V)
        g_resid = j_now - j_target

        normF = maximum(abs.(Fphys))
        rel_j = abs(g_resid) / j_target
        verbose && @info "galvanostatic Newton" iter normF rel_j V j_now

        if !isfinite(normF) || !isfinite(j_now)
            return (converged = false, iter = iter, normF = normF, u = u, V = V)
        end
        # Convergence: j matched to tolerance. Physics convergence is NOT required here —
        # the caller (_solve_galvanostatic!) will do a V-fixed physics refinement afterward.
        if rel_j < tol_rel_j
            return (converged = true, iter = iter, normF = normF, u = u, V = V)
        end

        # ∂F_phys/∂u — banded FD Jacobian
        banded_fd_jacobian!(Jsp, residual!, u, Fphys)
        @inbounds for i in 1:n; Jsp[i, i] += lambda_fix; end

        # Factor J (banded LU, reused for 2 solves)
        local LU_J
        try
            LU_J = lu(Jsp)
        catch e
            e isa LinearAlgebra.SingularException &&
                return (converged = false, iter = iter, normF = normF, u = u, V = V)
            rethrow(e)
        end

        # ∂F_phys/∂V — one extra residual eval
        build_residual(V + eps_V)(Fpert, u)
        dFdV = (Fpert .- Fphys) ./ eps_V

        # ∂j_total/∂u[6]  (log c_AN, electrode cell 1)
        eps6  = 1.0e-7 * max(abs(u[6]), 1.0)
        u_tmp = copy(u); u_tmp[6] += eps6
        dg_du6  = (j_total_fn(u_tmp, V) - j_now) / eps6

        # ∂j_total/∂u[10] (φ_l, electrode cell 1 — DOF 10 in v7)
        u_tmp = copy(u); u_tmp[10] += 1.0e-5
        dg_du10 = (j_total_fn(u_tmp, V) - j_now) / 1.0e-5

        # ∂j_total/∂V
        dg_dV = (j_total_fn(u, V + eps_V) - j_now) / eps_V

        # Bordered system solve via Schur complement:
        #   [J  c] [du]   [-Fphys]      c = dFdV
        #   [r' d] [dv] = [-g_resid]    r sparse (nonzeros at 6, 10), d = dg_dV
        #
        # p = J⁻¹ Fphys,  q = J⁻¹ c  (two banded LU solves)
        # dv = (-g_resid + r'p) / (d - r'q)
        # du = -(p + q·dv)
        p = LU_J \ Fphys
        q = LU_J \ dFdV

        r_dot_p = dg_du6 * p[6] + dg_du10 * p[10]
        r_dot_q = dg_du6 * q[6] + dg_du10 * q[10]
        S       = dg_dV - r_dot_q                    # Schur complement (scalar)

        abs(S) < 1.0e-20 &&
            return (converged = false, iter = iter, normF = normF, u = u, V = V)

        dv = (-g_resid + r_dot_p) / S
        du = -(p .+ q .* dv)

        # Step clamping (v7: 10 DOF/cell, φ_l at index 10)
        @inbounds for ix in 1:N_cells
            base = 10 * (ix - 1)
            for k in 1:9
                du[base + k] = clamp(du[base + k], -5.0, 5.0)
            end
            du[base + 10] = clamp(du[base + 10], -0.015, 0.015)
        end
        dv = clamp(dv, -0.10, 0.10)    # ≤ 100 mV per Newton step

        # When far from the current target (rel_j > 0.5), skip the line search and
        # accept the clamped step unconditionally. The line search cannot make useful
        # progress when j_stall << j_target because the fractional merit improvement
        # per step is too small to pass the strict-decrease test, yet any physics
        # disruption from a larger step pushes the combined merit above the starting
        # value. The 100 mV / 5-log / 15 mV step clamps already bound the step size.
        if rel_j > 0.5
            u_try = u .+ du
            V_try = V  + dv
            j_try = j_total_fn(u_try, V_try)
            if isfinite(j_try)
                u .= u_try
                V  = V_try
            end
        else
            # Physics-aware line search (used once j is within 50% of target):
            # When physics is already satisfied (normF < tol_phys), track only
            # |g|/j_target to avoid the deadlock where reducing |g| slightly
            # increases ||Fphys|| and the combined merit never decreases.
            merit_fn = (Fvec, g_val) ->
                maximum(abs.(Fvec)) < tol_phys ?
                abs(g_val) / j_target :
                sqrt(sum(abs2, Fvec) + (g_val / j_target)^2)

            normR2   = merit_fn(Fphys, g_resid)
            α        = 1.0
            accepted = false
            for _ in 0:max_backtrack
                u_try = u .+ α .* du
                V_try = V  + α  * dv
                build_residual(V_try)(Fpert, u_try)
                j_try      = j_total_fn(u_try, V_try)
                g_try      = j_try - j_target
                normR2_try = merit_fn(Fpert, g_try)
                if isfinite(normR2_try) && normR2_try < normR2
                    u .= u_try
                    V  = V_try
                    accepted = true
                    break
                end
                α *= 0.5
            end
            if !accepted
                verbose && @warn "galvanostatic: line search failed" iter normF
                return (converged = false, iter = iter, normF = normF, u = u, V = V)
            end
        end
    end

    return (converged = false, iter = max_iter,
            normF = maximum(abs.(Fphys)), u = u, V = V)
end

# ------------------------------------------------------------
# Newton continuation in V (§10.4)
# ------------------------------------------------------------
"""
    newton_continuation(u0, V_start, V_end, build_residual;
                        ds_init=0.05, ds_min=0.005, ds_max=0.20,
                        max_iter=25, verbose=false)
        -> Vector{Tuple{Float64, Vector{Float64}}}

Simple continuation sweeping V from V_start toward V_end (more negative).
`build_residual(V)` must return a closure `residual!(F, u)` for that V.
Adaptive: halves ds on failure, grows by 1.5× on fast (≤5 iter) success.

Step polarity: (V_end - V_start) < 0 ⇒ we step V by −ds each success.
"""
function newton_continuation(u0::Vector{Float64},
                             V_start::Float64, V_end::Float64,
                             build_residual::BF;
                             ds_init::Float64 = 0.05,
                             ds_min::Float64  = 1.0e-4,
                             ds_max::Float64  = 0.20,
                             max_iter::Int    = 60,
                             max_total_fail::Int = 200,
                             jacobian_mode::Symbol = :fd,
                             verbose::Bool    = false) where {BF<:Function}

    @assert V_end < V_start "Continuation expects V sweep toward more negative V"

    history = Tuple{Float64, Vector{Float64}}[]
    u  = copy(u0)
    V  = V_start
    ds = ds_init

    # always record starting point if we treat it as converged by caller
    # (caller should have already converged u0 at V_start)

    # First: converge at V_start (warm-start still wants true steady state).
    residual! = build_residual(V_start)
    res = newton_solve!(u, residual!;
                        max_iter = max_iter,
                        jacobian_mode = jacobian_mode,
                        verbose = verbose)
    if !res.converged
        error("Continuation: initial Newton at V = $V_start failed (|F| = $(res.normF)).")
    end
    push!(history, (V_start, copy(u)))
    verbose && @info "continuation: converged start" V=V_start iter=res.iter normF=res.normF

    # helper visible to the user so they can spell-check the struct
    total_failures = 0
    while V > V_end + 1e-12
        V_next = max(V - ds, V_end)
        residual! = build_residual(V_next)
        u_trial   = copy(u)
        res       = newton_solve!(u_trial, residual!;
                                  max_iter = max_iter,
                                  jacobian_mode = jacobian_mode,
                                  verbose = verbose)

        if res.converged
            verbose && @info "cont step OK" V=V_next ds=ds iter=res.iter normF=res.normF
            push!(history, (V_next, copy(u_trial)))
            u = u_trial
            V = V_next
            # Graduated adaptive growth/shrink based on Newton iter count.
            if res.iter ≤ 4
                ds = min(ds * 1.4, ds_max)
            elseif res.iter ≤ 10
                ds = min(ds * 1.1, ds_max)        # mild grow on moderate convergence
            elseif res.iter ≤ 20
                # stay (no-op)
            else
                ds = max(ds * 0.7, ds_min)        # shrink on slow convergence
            end
        else
            ds *= 0.3                              # aggressive shrink on failure
            total_failures += 1
            verbose && @info "cont step FAIL" V_target=V_next ds_new=ds normF=res.normF
            if ds < ds_min
                @warn "Continuation failure — ds below floor" V_last=V V_target=V_next ds ds_min
                break
            end
            if total_failures ≥ max_total_fail
                @warn "Continuation failure — too many total failures" V_last=V total_failures
                break
            end
        end
    end
    return history
end

# ============================================================
# Pseudo-arclength-style continuation parameterised by log₁₀(j_total)
# ============================================================
#
# Motivation:
#   Plain V-continuation stalls in the stiff regime past ~V = −2.3 V because
#   j_r(V) is nearly exponential — small V changes produce decades of current
#   change, and the continuation step size has to shrink below its floor to
#   keep Newton converging.  Stepping instead in log(j_total) auto-adapts:
#   uniform log-j increments correspond to ΔV ≈ RT/αF per decade.
#
# Augmented system (size n+1, where n = 9·N_cells):
#   [ F_physics(u, V) ]   = 0
#   [ log10(j_tot) − log10(j_tgt) ] = 0
#
# Unknowns: w = (u, V).  For the Jacobian we use the existing banded FD for
# ∂F/∂u, one extra residual eval for ∂F/∂V (dense column), and analytic/FD
# for the last row.  The last row has only 3 structurally nonzero entries:
# ∂g/∂u[6] (c_AN,1), ∂g/∂u[9] (φ_ℓ,1), and ∂g/∂V.
#
# ------------------------------------------------------------
"""
    newton_continuation_logj(u_start, V_start, V_end,
                             build_residual,
                             faradaic_current::Function;
                             n_steps_per_decade = 20,
                             max_iter = 40,
                             verbose = false)
        -> Vector{Tuple{Float64, Vector{Float64}}}

Continue the solution from V_start (already converged) to V_end (more negative)
in log-j steps.  `faradaic_current(u, V) -> j_total` computes the total current
density at the electrode from the current state and V.
"""
function newton_continuation_logj(u_start::Vector{Float64},
                                  V_start::Float64, V_end::Float64,
                                  build_residual::BF,
                                  faradaic_current::FJ;
                                  n_steps_per_decade::Int = 20,
                                  max_iter::Int = 40,
                                  tol::Float64  = 1.0e-8,
                                  verbose::Bool = false) where {BF<:Function, FJ<:Function}

    @assert V_end < V_start "log-j continuation expects V_end more negative than V_start"

    n = length(u_start)
    history = Tuple{Float64, Vector{Float64}}[(V_start, copy(u_start))]

    # initial j_total
    j_start = faradaic_current(u_start, V_start)
    @assert j_start > 0 "faradaic_current must be positive for log-j continuation"
    log10_j_start = log10(j_start)
    verbose && @info "logj cont: start" V=V_start j=j_start log10_j=log10_j_start

    # Pre-run a V-ladder estimate of j at V_end to know how many log-decades to cover.
    # For cathodic Tafel, j grows as V becomes MORE NEGATIVE, so d(log₁₀ j)/dV < 0.
    # Estimate: α_HER·F·log₁₀(e)/RT ≈ 0.4·38.94·0.4343 ≈ 6.8 (magnitude).
    dlogj_per_dV = -6.8
    log10_j_end  = log10_j_start + (V_end - V_start) * dlogj_per_dV
    total_decades = log10_j_end - log10_j_start
    N_steps = max(5, Int(round(abs(total_decades) * n_steps_per_decade)))
    verbose && @info "logj cont: target path" log10_j_start log10_j_end N_steps

    # Step schedule: uniform in log10(j).
    dlogj = (log10_j_end - log10_j_start) / N_steps

    u  = copy(u_start)
    V  = V_start
    Fphys = zeros(n)
    Fpert = zeros(n)

    for step in 1:N_steps
        log10_j_target = log10_j_start + step * dlogj

        # -------- augmented Newton loop --------
        converged_step = false
        for iter in 1:max_iter
            # --- physics residual at current (u, V) ---
            residual! = build_residual(V)
            residual!(Fphys, u)
            j_now = faradaic_current(u, V)
            g     = log10(max(j_now, 1e-300)) - log10_j_target

            R = vcat(Fphys, g)
            normR = maximum(abs.(R))
            verbose && iter == 1 && @info "logj step start" step log10_j_target V normR
            if normR < tol
                converged_step = true
                verbose && @info "logj step converged" step iter V normR
                break
            end

            # --- build augmented Jacobian ---
            #  top-left: banded FD (reuse existing)
            Jphys = build_banded_pattern(n, JAC_HALFBW)
            banded_fd_jacobian!(Jphys, (Fout, uu) -> begin
                                  residual!(Fout, uu)
                                  nothing
                              end, u, Fphys)
            # small LM damping (matches newton_solve!)
            @inbounds for i in 1:n
                Jphys[i, i] += 1.0e-10
            end

            # top-right column: ∂F_phys/∂V (FD, one residual eval)
            eps_V = 1e-5
            res_V_pert = build_residual(V + eps_V)
            res_V_pert(Fpert, u)
            dFdV = (Fpert .- Fphys) ./ eps_V

            # last row (∂g/∂u, ∂g/∂V):
            # g(u, V) = log10(j_tot) − log10_j_tgt
            # j_tot = Σ_r j_r depends on c_AN,surface (=u[6]) and φ_l,surface (=u[9]) and V.
            # Use FD.
            eps_u = 1e-7
            eps_phi = 1e-5
            # du[6] (c_AN log-conc, electrode cell 1)
            u_tmp = copy(u); u_tmp[6] += eps_u * max(abs(u[6]), 1.0)
            dg_du6 = (log10(max(faradaic_current(u_tmp, V), 1e-300)) - log10(j_now)) / (eps_u * max(abs(u[6]), 1.0))
            # du[10] (φ_l, electrode cell 1 — DOF 10 in v7, not 9)
            u_tmp = copy(u); u_tmp[10] += eps_phi
            dg_du9 = (log10(max(faradaic_current(u_tmp, V), 1e-300)) - log10(j_now)) / eps_phi
            # dV
            dg_dV = (log10(max(faradaic_current(u, V + eps_V), 1e-300)) - log10(j_now)) / eps_V

            # assemble augmented (n+1) × (n+1) sparse system
            # for simplicity, convert to dense — n+1 is ~901, dense LU is fast.
            Jaug = zeros(n + 1, n + 1)
            Jaug[1:n, 1:n] = Matrix(Jphys)
            Jaug[1:n, n+1] = dFdV
            Jaug[n+1, 6]   = dg_du6
            Jaug[n+1, 10]  = dg_du9   # φ_l is DOF 10 in v7
            Jaug[n+1, n+1] = dg_dV

            # solve
            dw = try
                -(Jaug \ R)
            catch e
                if e isa LinearAlgebra.SingularException
                    verbose && @warn "logj Newton: singular aug Jacobian" step iter
                    break
                else
                    rethrow(e)
                end
            end
            if !all(isfinite, dw)
                verbose && @warn "logj Newton: non-finite step" step iter
                break
            end

            # clamp per-DOF steps (v7: 10 DOF/cell, φ_l at DOF 10)
            N_cells = n ÷ 10
            @inbounds for ix in 1:N_cells
                base = 10 * (ix - 1)
                for k in 1:9
                    dw[base + k] = clamp(dw[base + k], -5.0, 5.0)
                end
                dw[base + 10] = clamp(dw[base + 10], -0.015, 0.015)
            end
            dw[n + 1] = clamp(dw[n + 1], -0.05, 0.05)   # ≤ 50 mV in V per Newton iter

            # L2 backtracking line search
            normR2 = sqrt(sum(abs2, R))
            α = 1.0
            accepted = false
            for _ in 0:10
                u_try = u .+ α .* dw[1:n]
                V_try = V + α * dw[n + 1]
                residual_try = build_residual(V_try)
                residual_try(Fpert, u_try)
                g_try = log10(max(faradaic_current(u_try, V_try), 1e-300)) - log10_j_target
                R_try = vcat(Fpert, g_try)
                normR2_try = sqrt(sum(abs2, R_try))
                if isfinite(normR2_try) && normR2_try < normR2
                    u .= u_try
                    V  = V_try
                    accepted = true
                    break
                end
                α *= 0.5
            end
            if !accepted
                verbose && @warn "logj Newton: line search failed" step iter
                break
            end
        end

        if converged_step
            push!(history, (V, copy(u)))
        else
            @warn "logj continuation: step failed" step log10_j_target V
            break
        end

        V <= V_end + 1e-12 && break
    end

    return history
end

end # module
