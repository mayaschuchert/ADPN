# run_stage4_seq.jl — Sequential decoupled kinetics fit (Bloomquist SI §S1.11)
#
# Fits one reaction at a time while holding the other reactions' kinetic
# parameters (j0, alpha_c, n) fixed at their current best estimates.
# Each sub-problem is 3 parameters: (log10 j0_r, alpha_c_r, n_r).
# solve_at_j runs normally (no FIXED_J); FE_r comes from the converged PDE.
#
# Fit order: ADN → PN → TCH.  HER frozen at v6.x values throughout.
# Passes: up to MAX_PASSES Gauss-Seidel cycles or until all channels converge.
# Output: output/stage4_seq/
using Printf
using Dates
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics

const DATA_FILE    = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const OUT_DIR      = joinpath(@__DIR__, "output", "stage4_seq")
const OUT_DATA_DIR = joinpath(OUT_DIR, "data")
const OUT_LOG_DIR  = joinpath(OUT_DIR, "logs")
isdir(OUT_DATA_DIR) || mkpath(OUT_DATA_DIR)
isdir(OUT_LOG_DIR)  || mkpath(OUT_LOG_DIR)

const MAX_PASSES = 3

# ---------- CSV loader (same as run_stage4v3.jl) ----------
function load_bloomquist(path::String)
    raw, hdr = readdlm(path, ','; header = true)
    hdr = vec(hdr)
    col(name) = findfirst(==(string(name)), hdr)
    rows = BloomquistRow[]
    for i in 1:size(raw, 1)
        push!(rows, BloomquistRow(
            string(raw[i, col("table")]),
            Float64(raw[i, col("gap_mm")]),
            Float64(raw[i, col("Q_total_mL_min")]),
            Float64(raw[i, col("j_mA_cm2")]),
            Float64(raw[i, col("phi_AN")]),
            Float64(raw[i, col("Q_aq_mL_min")]),
            Float64(raw[i, col("Q_org_mL_min")]),
            Float64(raw[i, col("We_aq")]),
            Float64(raw[i, col("We_org")]),
            Float64(raw[i, col("FE_ADN_pct")]),
            Float64(raw[i, col("FE_TCH_pct")]),
            Float64(raw[i, col("FE_PN_pct")]),
            Float64(raw[i, col("PR_ADN_kg_cm2_h")]),
            Float64(raw[i, col("EP_ADN_kg_kWh")]),
            NaN, NaN, NaN, NaN
        ))
    end
    return rows
end

# ---------- Single-reaction 1-channel residual ----------
# reaction_id: 1=ADN, 2=PN, 4=TCH  (3=HER frozen, never fit)
# theta3 = [log10(j0_r), alpha_c_r, n_r]
# j0_cur/ac_cur/n_cur = current best full-parameter tuples for ALL reactions.
# sol_cache: maps row key → (V_cathode, state) from the most recent successful
#   solve. Passed as V_hint/u_hint so solve_at_j skips the walk for FD columns.
function seq_residuals!(F::Vector{Float64}, theta3::Vector{Float64},
                        ctx::FitContext, sel::Vector{Int}, reaction_id::Int,
                        j0_cur::NTuple{4,Float64},
                        ac_cur::NTuple{4,Float64},
                        n_cur::NTuple{3,Float64};
                        sol_cache::Dict{Any,Tuple{Float64,Vector{Float64}}} =
                            Dict{Any,Tuple{Float64,Vector{Float64}}}())
    nfail = 0
    for (k, idx) in pairs(sel)
        r = ctx.rows[idx]

        # Patch theta3 into the full kinetics tuples for reaction_id only
        j0_ov = j0_cur; ac_ov = ac_cur; n_ov = n_cur
        if reaction_id == 1
            j0_ov = (10.0^theta3[1], j0_ov[2], j0_ov[3], j0_ov[4])
            ac_ov = (theta3[2],      ac_ov[2], ac_ov[3], ac_ov[4])
            n_ov  = (theta3[3],      n_ov[2],  n_ov[3])
        elseif reaction_id == 2
            j0_ov = (j0_ov[1], 10.0^theta3[1], j0_ov[3], j0_ov[4])
            ac_ov = (ac_ov[1], theta3[2],       ac_ov[3], ac_ov[4])
            n_ov  = (n_ov[1],  theta3[3],       n_ov[3])
        else  # reaction_id == 4 (TCH)
            j0_ov = (j0_ov[1], j0_ov[2], j0_ov[3], 10.0^theta3[1])
            ac_ov = (ac_ov[1], ac_ov[2], ac_ov[3], theta3[2])
            n_ov  = (n_ov[1],  n_ov[2],  theta3[3])
        end

        key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
        δkey = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[δkey]
        u_w  = get(ctx.warm_by_key, key,
                   make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN))

        # Pull (V_cathode, state) hint if available
        V_h, u_h = if haskey(sol_cache, key)
            sol_cache[key]
        else
            (NaN, nothing)
        end

        result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                            mesh, u_w, ctx.c_eq;
                            j0 = j0_ov, alpha_c = ac_ov, n_orders = n_ov,
                            V_hint = V_h, u_hint = u_h,
                            ds_min = 1.0e-5, newton_max_iter = 120)

        # Fallback: if walk V_hi was below V_target (j(V_hi) ≥ j_target), redo without hint
        if !result.converged && occursin("already ≥ j_target", result.note)
            result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                                mesh, u_w, ctx.c_eq;
                                j0 = j0_ov, alpha_c = ac_ov, n_orders = n_ov,
                                ds_min = 1.0e-5, newton_max_iter = 120)
        end

        if result.converged
            sol_cache[key] = (result.V_cathode, result.state)
            ctx.warm_by_key[key] = result.state
            FE_mod = (reaction_id == 1 ? result.FE_ADN_pct :
                      reaction_id == 2 ? result.FE_PN_pct  : result.FE_TCH_pct)
            FE_obs = (reaction_id == 1 ? r.FE_ADN_pct :
                      reaction_id == 2 ? r.FE_PN_pct  : r.FE_TCH_pct)
            F[k] = FE_mod - FE_obs
        else
            if nfail == 0
                @printf("[seq_residuals! FAIL] row=%d  note: %s\n", idx, result.note)
                flush(stdout)
            end
            F[k] = 1.0e3
            nfail += 1
        end
    end
    return nfail
end

# ---------- Minimal 3-parameter LM ----------
function lm3(f!, theta0::Vector{Float64},
             lb::Vector{Float64}, ub::Vector{Float64};
             nrows::Int,
             lambda0::Float64 = 1e-2,
             max_iter::Int    = 60,
             tol_rel::Float64 = 1e-2,
             fd_step::Float64 = 1e-4,
             reaction_name::String = "?",
             verbose::Bool = true)

    theta = clamp.(copy(theta0), lb, ub)
    F     = zeros(nrows)
    nfail = f!(F, theta)
    loss  = dot(F, F)
    lambda = lambda0

    if verbose
        @printf("[%s] iter  0  loss=%.4e  λ=%.2e  nfail=%d\n",
                reaction_name, loss, lambda, nfail)
        flush(stdout)
    end

    J      = zeros(nrows, 3)
    F_pert = zeros(nrows)
    best_theta = copy(theta)
    best_loss  = loss

    for it in 1:max_iter
        for k in 1:3
            dθ = max(fd_step, fd_step * abs(theta[k]))
            θk = theta[k]
            theta[k] = clamp(θk + dθ, lb[k], ub[k])
            f!(F_pert, theta)
            J[:, k] = (F_pert .- F) ./ (theta[k] - θk + eps())
            theta[k] = θk
        end

        JtJ = J' * J
        H   = JtJ + lambda * Diagonal(diag(JtJ) .+ 1e-12)
        Δ   = -(H \ (J' * F))

        theta_trial = clamp.(theta .+ Δ, lb, ub)
        F_trial = zeros(nrows)
        nfail_t = f!(F_trial, theta_trial)
        loss_t  = dot(F_trial, F_trial)

        if loss_t < loss
            rel = (loss - loss_t) / max(loss, 1e-12)
            theta  = theta_trial
            F     .= F_trial
            loss   = loss_t
            nfail  = nfail_t
            lambda = max(lambda * 0.5, 1e-8)
            if loss < best_loss; best_theta = copy(theta); best_loss = loss; end
            if verbose
                @printf("[%s] iter %2d  loss=%.4e  λ=%.2e  nfail=%d  ACCEPT Δrel=%.3e\n",
                        reaction_name, it, loss, lambda, nfail, rel)
                @printf("         theta=[%.4f, %.4f, %.4f]  → j0=%.3e  ac=%.3f  n=%.3f\n",
                        theta[1], theta[2], theta[3],
                        10.0^theta[1], theta[2], theta[3])
                flush(stdout)
            end
            if rel < tol_rel
                return (theta = best_theta, loss = best_loss, iter = it,
                        converged = true, note = "rel drop < tol_rel")
            end
        else
            lambda = min(lambda * 4.0, 1e4)
            if verbose
                @printf("[%s] iter %2d  loss=%.4e  λ=%.2e  nfail=%d  REJECT\n",
                        reaction_name, it, loss, lambda, nfail)
                flush(stdout)
            end
            if lambda > 1e3
                return (theta = best_theta, loss = best_loss, iter = it,
                        converged = false, note = "lambda stuck")
            end
        end
    end
    return (theta = best_theta, loss = best_loss, iter = max_iter,
            converged = false, note = "max_iter")
end

function rmse1(F::Vector{Float64})
    return sqrt(dot(F, F) / max(length(F), 1))
end

function _pin(val, lb, ub; tol=1e-3)
    span = ub - lb
    (val - lb) ≤ tol * span && return "@LB"
    (ub - val) ≤ tol * span && return "@UB"
    return "   "
end

# ---------- Main ----------
function main()
    println("=" ^ 72)
    println(" Stage 4-seq — Sequential decoupled fit (Bloomquist §S1.11)")
    println(" $(now())")
    println("=" ^ 72)
    flush(stdout)

    rows_raw = load_bloomquist(DATA_FILE)
    @printf("Loaded %d rows\n", length(rows_raw))

    sel_core     = select_core(rows_raw)
    sel_extended = select_extended(rows_raw)
    sel_holdout  = select_holdout(rows_raw)
    @printf("Core=%d  Extended=%d  Holdout=%d\n\n",
            length(sel_core), length(sel_extended), length(sel_holdout))
    flush(stdout)

    ctx = build_context(rows_raw, sel_core)
    @printf("Context built: %d meshes, %d warm starts\n\n",
            length(ctx.mesh_by_delta), length(ctx.warm_by_key))
    flush(stdout)

    # ---- Parameter bounds (per reaction, 3 params each) ----
    # theta3 = [log10(j0_r), alpha_c_r, n_r]
    LB_ADN = [-6.0, 0.30, 0.5];  UB_ADN = [-1.0, 0.70, 3.0]
    LB_PN  = [-6.0, 0.30, 0.5];  UB_PN  = [-1.0, 0.70, 2.0]
    LB_TCH = [-6.0, 0.30, 1.0];  UB_TCH = [-1.0, 0.70, 3.0]

    # ---- Initial guess from THETA0 ----
    # THETA0 layout: [log10j0_ADN, log10j0_PN, log10j0_TCH,
    #                  ac_ADN, ac_PN, ac_TCH, n_ADN, n_PN, n_TCH]
    t0_ADN = [THETA0[1], THETA0[4], THETA0[7]]
    t0_PN  = [THETA0[2], THETA0[5], THETA0[8]]
    t0_TCH = [THETA0[3], THETA0[6], THETA0[9]]

    # Current best full tuples (updated after each reaction fit)
    j0_cur = (10.0^THETA0[1], 10.0^THETA0[2], J0_3_FROZEN, 10.0^THETA0[3])
    ac_cur = (THETA0[4], THETA0[5], ALPHA_C3_FROZEN, THETA0[6])
    n_cur  = (THETA0[7], THETA0[8], THETA0[9])

    best_ADN = copy(t0_ADN)
    best_PN  = copy(t0_PN)
    best_TCH = copy(t0_TCH)

    nrows = length(sel_core)

    # ---- Pre-solve at THETA0: populate sol_cache and warm starts ----
    # This one-time walk from V_hi=-0.8 makes every subsequent FD eval fast:
    # FD columns use V_hint/u_hint from sol_cache → 1-2 Newton steps instead of
    # the full 1.1 V walk, giving ~10-50x speedup per residual evaluation.
    println("\n[diag] Pre-solving core rows at THETA0 (populating sol_cache)...")
    flush(stdout)
    sol_cache = Dict{Any,Tuple{Float64,Vector{Float64}}}()
    j0_init = (10.0^THETA0[1], 10.0^THETA0[2], J0_3_FROZEN, 10.0^THETA0[3])
    ac_init = (THETA0[4], THETA0[5], ALPHA_C3_FROZEN, THETA0[6])
    n_init  = (THETA0[7], THETA0[8], THETA0[9])
    nfail_init = 0
    loss_init  = 0.0
    for (_, idx) in pairs(sel_core)
        r    = ctx.rows[idx]
        key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
        δkey = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[δkey]
        u_w  = get(ctx.warm_by_key, key,
                   make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN))
        result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                            mesh, u_w, ctx.c_eq;
                            j0 = j0_init, alpha_c = ac_init, n_orders = n_init,
                            ds_min = 1.0e-5, newton_max_iter = 120)
        if result.converged
            sol_cache[key] = (result.V_cathode, result.state)
            ctx.warm_by_key[key] = result.state
            loss_init += (result.FE_ADN_pct - r.FE_ADN_pct)^2 +
                         (result.FE_PN_pct  - r.FE_PN_pct)^2  +
                         (result.FE_TCH_pct - r.FE_TCH_pct)^2
        else
            nfail_init += 1
            loss_init  += 3.0 * 1.0e6
        end
    end
    @printf("[diag] THETA0 pre-solve: nfail=%d / %d  loss=%.4e  sol_cache=%d entries\n\n",
            nfail_init, nrows, loss_init, length(sol_cache))
    flush(stdout)

    for pass in 1:MAX_PASSES
        println("\n" * "─" ^ 72)
        @printf(" Gauss-Seidel pass %d / %d\n", pass, MAX_PASSES)
        println("─" ^ 72)
        flush(stdout)

        # ---- Fit ADN ----
        println("\n[ADN] Fixing PN/HER/TCH to experimental; fitting j0_ADN, α_ADN, n_ADN")
        flush(stdout)
        f_ADN! = (F, t) -> seq_residuals!(F, t, ctx, sel_core, 1,
                                           j0_cur, ac_cur, n_cur;
                                           sol_cache = sol_cache)
        res_ADN = lm3(f_ADN!, best_ADN, LB_ADN, UB_ADN;
                      nrows = nrows, reaction_name = "ADN")
        best_ADN = res_ADN.theta
        j0_cur = (10.0^best_ADN[1], j0_cur[2], j0_cur[3], j0_cur[4])
        ac_cur = (best_ADN[2], ac_cur[2], ac_cur[3], ac_cur[4])
        n_cur  = (best_ADN[3], n_cur[2],  n_cur[3])
        @printf("[ADN] done: j0=%.3e %s  α=%.3f %s  n=%.3f %s  RMSE=%.2f pp\n",
                j0_cur[1], _pin(best_ADN[1], LB_ADN[1], UB_ADN[1]),
                ac_cur[1], _pin(best_ADN[2], LB_ADN[2], UB_ADN[2]),
                n_cur[1],  _pin(best_ADN[3], LB_ADN[3], UB_ADN[3]),
                sqrt(res_ADN.loss / nrows))
        flush(stdout)

        # ---- Fit PN ----
        println("\n[PN]  Fixing ADN(fitted)/HER/TCH to experimental; fitting j0_PN, α_PN, n_PN")
        flush(stdout)
        f_PN! = (F, t) -> seq_residuals!(F, t, ctx, sel_core, 2,
                                          j0_cur, ac_cur, n_cur;
                                          sol_cache = sol_cache)
        res_PN = lm3(f_PN!, best_PN, LB_PN, UB_PN;
                     nrows = nrows, reaction_name = "PN")
        best_PN = res_PN.theta
        j0_cur = (j0_cur[1], 10.0^best_PN[1], j0_cur[3], j0_cur[4])
        ac_cur = (ac_cur[1], best_PN[2], ac_cur[3], ac_cur[4])
        n_cur  = (n_cur[1],  best_PN[3], n_cur[3])
        @printf("[PN]  done: j0=%.3e %s  α=%.3f %s  n=%.3f %s  RMSE=%.2f pp\n",
                j0_cur[2], _pin(best_PN[1], LB_PN[1], UB_PN[1]),
                ac_cur[2], _pin(best_PN[2], LB_PN[2], UB_PN[2]),
                n_cur[2],  _pin(best_PN[3], LB_PN[3], UB_PN[3]),
                sqrt(res_PN.loss / nrows))
        flush(stdout)

        # ---- Fit TCH ----
        println("\n[TCH] Fixing ADN(fitted)/PN(fitted)/HER to experimental; fitting j0_TCH, α_TCH, n_TCH")
        flush(stdout)
        f_TCH! = (F, t) -> seq_residuals!(F, t, ctx, sel_core, 4,
                                           j0_cur, ac_cur, n_cur;
                                           sol_cache = sol_cache)
        res_TCH = lm3(f_TCH!, best_TCH, LB_TCH, UB_TCH;
                      nrows = nrows, reaction_name = "TCH")
        best_TCH = res_TCH.theta
        j0_cur = (j0_cur[1], j0_cur[2], j0_cur[3], 10.0^best_TCH[1])
        ac_cur = (ac_cur[1], ac_cur[2], ac_cur[3], best_TCH[2])
        n_cur  = (n_cur[1],  n_cur[2],  best_TCH[3])
        @printf("[TCH] done: j0=%.3e %s  α=%.3f %s  n=%.3f %s  RMSE=%.2f pp\n",
                j0_cur[4], _pin(best_TCH[1], LB_TCH[1], UB_TCH[1]),
                ac_cur[4], _pin(best_TCH[2], LB_TCH[2], UB_TCH[2]),
                n_cur[3],  _pin(best_TCH[3], LB_TCH[3], UB_TCH[3]),
                sqrt(res_TCH.loss / nrows))
        flush(stdout)
    end

    # ---- Forward apply to Extended + Holdout ----
    println("\n" * "=" ^ 72)
    println(" Forward apply to Extended + Holdout")
    println("=" ^ 72)
    flush(stdout)

    function eval_set(sel, label)
        ctx_s = build_context(rows_raw, sel; warm_init = copy(ctx.warm_by_key))
        F_adn = zeros(length(sel)); F_pn = zeros(length(sel)); F_tch = zeros(length(sel))
        nfail = 0
        for (k, idx) in pairs(sel)
            r = ctx_s.rows[idx]
            key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
            δkey = round(r.delta_lev_m * 1e6; digits = 1)
            mesh = ctx_s.mesh_by_delta[δkey]
            u_w  = get(ctx_s.warm_by_key, key,
                       make_initial_guess(ctx_s.N_mesh, ctx_s.c_eq, r.phi_AN))
            res = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                             mesh, u_w, ctx_s.c_eq;
                             j0 = j0_cur, alpha_c = ac_cur, n_orders = n_cur)
            if res.converged
                ctx_s.warm_by_key[key] = res.state
                F_adn[k] = res.FE_ADN_pct - r.FE_ADN_pct
                F_pn[k]  = res.FE_PN_pct  - r.FE_PN_pct
                F_tch[k] = res.FE_TCH_pct - r.FE_TCH_pct
            else
                nfail += 1
            end
        end
        @printf("  %-10s  ADN=%.2f pp  PN=%.2f pp  TCH=%.2f pp  nfail=%d\n",
                label, rmse1(F_adn), rmse1(F_pn), rmse1(F_tch), nfail)
        flush(stdout)
        return (F_adn, F_pn, F_tch)
    end

    F_core_adn = zeros(length(sel_core)); F_core_pn = zeros(length(sel_core))
    F_core_tch = zeros(length(sel_core))
    for (k, idx) in pairs(sel_core)
        r = ctx.rows[idx]
        key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
        δkey = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[δkey]
        u_w  = get(ctx.warm_by_key, key,
                   make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN))
        res = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                         mesh, u_w, ctx.c_eq;
                         j0 = j0_cur, alpha_c = ac_cur, n_orders = n_cur)
        if res.converged
            ctx.warm_by_key[key] = res.state
            F_core_adn[k] = res.FE_ADN_pct - r.FE_ADN_pct
            F_core_pn[k]  = res.FE_PN_pct  - r.FE_PN_pct
            F_core_tch[k] = res.FE_TCH_pct - r.FE_TCH_pct
        end
    end

    println("\nFinal RMSE (joint model eval — no fixed-j):")
    @printf("  %-10s  ADN=%.2f pp  PN=%.2f pp  TCH=%.2f pp\n",
            "Core", rmse1(F_core_adn), rmse1(F_core_pn), rmse1(F_core_tch))
    eval_set(sel_extended, "Extended")
    eval_set(sel_holdout,  "Holdout")

    # ---- Decision gates ----
    println("\n--- Decision gates ---")
    flush(stdout)
    function gate(name, val, thr)
        ok = val < thr
        @printf("  [%s] %-42s  %.2f pp  (< %.0f pp)\n",
                ok ? "PASS" : "FAIL", name, val, thr)
        flush(stdout)
        return ok
    end
    g1 = gate("Core FE_ADN RMSE",     rmse1(F_core_adn), 8.0)
    g2 = gate("Core FE_PN RMSE",      rmse1(F_core_pn),  5.0)
    g3 = gate("Core FE_TCH RMSE",     rmse1(F_core_tch), 4.0)
    println("\n  All gates passed: ", all((g1, g2, g3)))

    # ---- Save results ----
    open(joinpath(OUT_DATA_DIR, "stage4_seq_fitted_theta.txt"), "w") do io
        @printf(io, "# Sequential decoupled fit — %s\n", now())
        @printf(io, "j0_ADN      = %.6e\n", j0_cur[1])
        @printf(io, "j0_PN       = %.6e\n", j0_cur[2])
        @printf(io, "j0_HER      = %.6e  # frozen\n", J0_3_FROZEN)
        @printf(io, "j0_TCH      = %.6e\n", j0_cur[4])
        @printf(io, "alpha_c_ADN = %.6f\n", ac_cur[1])
        @printf(io, "alpha_c_PN  = %.6f\n", ac_cur[2])
        @printf(io, "alpha_c_HER = %.6f  # frozen\n", ALPHA_C3_FROZEN)
        @printf(io, "alpha_c_TCH = %.6f\n", ac_cur[4])
        @printf(io, "n_ADN       = %.6f\n", n_cur[1])
        @printf(io, "n_PN        = %.6f\n", n_cur[2])
        @printf(io, "n_TCH       = %.6f\n", n_cur[3])
    end

    println("\n" * "=" ^ 72)
    println(" Done.  Results in $(OUT_DATA_DIR)")
    println("=" ^ 72)
    flush(stdout)
end

main()
