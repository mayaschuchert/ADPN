# -----------------------------------------------------------------------------
# run_stage4v4.jl — Bloomquist kinetics fit, v7 §22: Lees S1.11 sequential
# per-branch fit on top of the new galvanostatic continuation solver.
#
# Changes from run_stage4v3.jl (joint 9-param LM, HER frozen):
#   - 4 sub-fits in order: HER → ADN → PN → TCH
#   - All non-active branches clamped to experimental partial currents via
#     Kinetics.FIXED_J (S1.11 prescription)
#   - HER is a fitted parameter (j0_HER, αc_HER), no longer frozen
#   - Inner solver: GalvContinuation.solve_at_j_cont (anchor + log-j sweep,
#     bordered (n+1)×(n+1) galvanostatic Newton at each step)
#   - Loss per sub-fit: Σ (log10(j_branch_model + 1) − log10(j_branch_obs + 1))²
#   - Optional outer iteration (default 2) to retighten cross-coupling
#   - δ_lev still NOT a fitted parameter — frozen at Lévêque value per (gap, Q)
#
# Pipeline:
#   Stage 4a: sequential fit on Core
#   Stage 4b: forward apply converged θ to Extended + Holdout (no re-fit) via
#             solve_at_j_cont, no clamping
#   Output  : output/stage4v4/
# -----------------------------------------------------------------------------
using Printf
using Dates
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.FitKineticsSeq
using .ADPN_EHD.GalvContinuation
using .ADPN_EHD.Kinetics

# ---------- Paths ----------
const DATA_FILE    = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const OUT_DIR      = joinpath(@__DIR__, "output", "stage4v4")
const OUT_DATA_DIR = joinpath(OUT_DIR, "data")
const OUT_LOG_DIR  = joinpath(OUT_DIR, "logs")
isdir(OUT_DATA_DIR) || mkpath(OUT_DATA_DIR)
isdir(OUT_LOG_DIR)  || mkpath(OUT_LOG_DIR)

# ---------- CSV loader (same as v3) ----------
function load_bloomquist(path::String)
    raw, hdr = readdlm(path, ','; header = true)
    hdr = vec(hdr)
    col(name) = findfirst(==(string(name)), hdr)
    @assert col("table") !== nothing "missing column: table"

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

# ---------- 4-channel residual writer (FE_ADN, FE_PN, FE_HER, FE_TCH) ----------
function write_residuals(path::String, rows::Vector{BloomquistRow},
                          sel::Vector{Int}, results::Vector)
    open(path, "w") do io
        println(io, "table,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN," *
                    "FE_ADN_obs,FE_ADN_model,FE_ADN_resid_pp," *
                    "FE_PN_obs,FE_PN_model,FE_PN_resid_pp," *
                    "FE_HER_obs,FE_HER_model,FE_HER_resid_pp," *
                    "FE_TCH_obs,FE_TCH_model,FE_TCH_resid_pp," *
                    "V_cathode,j_total_model,note")
        for (n, idx) in pairs(sel)
            r = rows[idx]
            res = results[n]
            fe_her_obs = max(0.0, 100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct))
            if res.converged
                @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.2f,%s\n",
                        r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.FE_ADN_pct, res.FE_ADN_pct, res.FE_ADN_pct - r.FE_ADN_pct,
                        r.FE_PN_pct,  res.FE_PN_pct,  res.FE_PN_pct  - r.FE_PN_pct,
                        fe_her_obs,   res.FE_HER_pct, res.FE_HER_pct - fe_her_obs,
                        r.FE_TCH_pct, res.FE_TCH_pct, res.FE_TCH_pct - r.FE_TCH_pct,
                        res.V_cathode, res.j_total / 10.0, res.note)
            else
                @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,NA,NA,%.2f,NA,NA,%.2f,NA,NA,%.2f,NA,NA,NA,NA,FAIL: %s\n",
                        r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.FE_ADN_pct,
                        r.FE_PN_pct,
                        fe_her_obs,
                        r.FE_TCH_pct,
                        res.note)
            end
        end
    end
end

function rmse_pp(rows::Vector{BloomquistRow}, sel::Vector{Int},
                  results::Vector, channel::Symbol)
    s = 0.0; n = 0
    for (i, idx) in pairs(sel)
        r = rows[idx]
        res = results[i]
        res.converged || continue
        if channel === :ADN
            s += (res.FE_ADN_pct - r.FE_ADN_pct)^2
        elseif channel === :PN
            s += (res.FE_PN_pct - r.FE_PN_pct)^2
        elseif channel === :TCH
            s += (res.FE_TCH_pct - r.FE_TCH_pct)^2
        elseif channel === :HER
            fe_her_obs = max(0.0, 100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct))
            s += (res.FE_HER_pct - fe_her_obs)^2
        else
            error("bad channel: $channel")
        end
        n += 1
    end
    return n == 0 ? NaN : sqrt(s / n)
end

# ---------- Forward apply (no clamping) ----------
# Reads derived fields (delta_lev_m, j_target_A_m2) from ctx.rows.
function forward_apply(sel::Vector{Int}, ctx, st::SeqFitState;
                        verbose::Bool = false)
    j0_use = (st.j0_ADN, st.j0_PN, st.j0_HER, st.j0_TCH)
    ac_use = (st.ac_ADN, st.ac_PN, st.ac_HER, st.ac_TCH)
    n_use  = (st.n_ADN, st.n_PN, st.n_TCH)

    results = Vector{ContResult}(undef, length(sel))
    nfail = 0
    for (i, idx) in pairs(sel)
        # Use ctx.rows (built / has delta_lev_m populated) instead of input rows
        r = ctx.rows[idx]
        key = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[key]
        result = solve_at_j_cont(r.j_target_A_m2, r.phi_AN, mesh, ctx.c_eq;
                                  j0       = j0_use,
                                  alpha_c  = ac_use,
                                  n_orders = n_use,
                                  fixed_j  = nothing,    # no clamping in forward apply
                                  warm_u   = nothing,
                                  verbose  = false)
        results[i] = result
        if !result.converged
            nfail += 1
            verbose && @warn "forward FAIL" row=idx note=result.note
        end
    end
    return results, nfail
end

# ---------- Pretty-print helpers ----------
function _print_state(st::SeqFitState)
    @printf("  j0_ADN = %.3e A/m²    ac_ADN = %.3f    n_ADN = %.3f\n",
            st.j0_ADN, st.ac_ADN, st.n_ADN)
    @printf("  j0_PN  = %.3e A/m²    ac_PN  = %.3f    n_PN  = %.3f\n",
            st.j0_PN,  st.ac_PN,  st.n_PN)
    @printf("  j0_HER = %.3e A/m²    ac_HER = %.3f    (fitted, was frozen in v3)\n",
            st.j0_HER, st.ac_HER)
    @printf("  j0_TCH = %.3e A/m²    ac_TCH = %.3f    n_TCH = %.3f\n",
            st.j0_TCH, st.ac_TCH, st.n_TCH)
end

# ---------- Main ----------
function main()
    println("=" ^ 72)
    println(" Stage 4v4 — Lees S1.11 sequential fit + galv-continuation solver")
    println(" $(now())")
    println("=" ^ 72)

    rows_raw = load_bloomquist(DATA_FILE)
    @printf("Loaded %d rows from %s\n", length(rows_raw), DATA_FILE)

    sel_core_idx     = select_core(rows_raw)
    sel_extended_idx = select_extended(rows_raw)
    sel_holdout_idx  = select_holdout(rows_raw)
    @printf("Core      : %3d rows\n", length(sel_core_idx))
    @printf("Extended  : %3d rows\n", length(sel_extended_idx))
    @printf("Holdout   : %3d rows\n", length(sel_holdout_idx))

    # ---- Stage 4a: sequential fit on Core ----
    println("\n--- Stage 4a: Sequential S1.11 fit on Core ---")
    ctx_core = build_context(rows_raw, sel_core_idx)
    @printf("Pre-cached %d unique-δ meshes\n", length(ctx_core.mesh_by_delta))

    # Warm-start from v3 converged params (not cold seed) — improvement #2.
    # Each per-branch fit then refines instead of searching from scratch.
    st, history = run_sequential_fit(ctx_core;
                                     warm_start      = SeqFitState_from_v3(),
                                     max_outer       = 5,
                                     branch_max_iter = 30,
                                     lambda_stuck    = 1.0e5,
                                     verbose         = true)

    println("\nFitted parameters (after sequential fit):")
    _print_state(st)

    # Build results vector at the converged state to compute RMSEs/residuals
    println("\n--- Stage 4a: re-evaluating Core at converged state (no clamping) ---")
    res_core, nf_core = forward_apply(sel_core_idx, ctx_core, st)
    rmse_core_adn = rmse_pp(rows_raw, sel_core_idx, res_core, :ADN)
    rmse_core_pn  = rmse_pp(rows_raw, sel_core_idx, res_core, :PN)
    rmse_core_her = rmse_pp(rows_raw, sel_core_idx, res_core, :HER)
    rmse_core_tch = rmse_pp(rows_raw, sel_core_idx, res_core, :TCH)
    @printf("Core RMSE (no clamping) — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_HER: %.2f pp,  FE_TCH: %.2f pp,  nfail=%d\n",
            rmse_core_adn, rmse_core_pn, rmse_core_her, rmse_core_tch, nf_core)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4a_core_residuals.csv"),
                    rows_raw, sel_core_idx, res_core)

    # Save fitted-θ
    open(joinpath(OUT_DATA_DIR, "stage4v4_fitted_theta.txt"), "w") do io
        @printf(io, "# v7 Stage 4v4 — Lees S1.11 sequential fit (HER not frozen, TCH active)\n")
        @printf(io, "# date: %s\n", now())
        @printf(io, "# Core rows: %d\n", length(sel_core_idx))
        @printf(io, "# RMSE FE_ADN: %.4f pp\n", rmse_core_adn)
        @printf(io, "# RMSE FE_PN:  %.4f pp\n", rmse_core_pn)
        @printf(io, "# RMSE FE_HER: %.4f pp (HER fitted, not frozen)\n", rmse_core_her)
        @printf(io, "# RMSE FE_TCH: %.4f pp\n", rmse_core_tch)
        @printf(io, "j0_ADN   = %.6e\n", st.j0_ADN)
        @printf(io, "j0_PN    = %.6e\n", st.j0_PN)
        @printf(io, "j0_HER   = %.6e\n", st.j0_HER)
        @printf(io, "j0_TCH   = %.6e\n", st.j0_TCH)
        @printf(io, "alpha_c_ADN = %.6f\n", st.ac_ADN)
        @printf(io, "alpha_c_PN  = %.6f\n", st.ac_PN)
        @printf(io, "alpha_c_HER = %.6f\n", st.ac_HER)
        @printf(io, "alpha_c_TCH = %.6f\n", st.ac_TCH)
        @printf(io, "n_ADN = %.6f\n", st.n_ADN)
        @printf(io, "n_PN  = %.6f\n", st.n_PN)
        @printf(io, "n_TCH = %.6f\n", st.n_TCH)
    end

    # Save sub-fit history
    open(joinpath(OUT_DATA_DIR, "stage4v4_seq_history.csv"), "w") do io
        println(io, "outer,branch,loss_log10j2,iter,note")
        for h in history
            @printf(io, "%d,%s,%.6e,%d,%s\n", h.outer, h.branch,
                    h.loss, h.iter, h.note)
        end
    end

    # ---- Stage 4b: forward apply ----
    println("\n--- Stage 4b: forward apply θ to Extended + Holdout ---")

    ctx_ext = build_context(rows_raw, sel_extended_idx;
                            warm_init = copy(ctx_core.warm_by_key))
    res_ext, nf_ext = forward_apply(sel_extended_idx, ctx_ext, st)
    rmse_ext_adn = rmse_pp(rows_raw, sel_extended_idx, res_ext, :ADN)
    rmse_ext_pn  = rmse_pp(rows_raw, sel_extended_idx, res_ext, :PN)
    rmse_ext_tch = rmse_pp(rows_raw, sel_extended_idx, res_ext, :TCH)
    @printf("Extended RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp,  nfail=%d\n",
            rmse_ext_adn, rmse_ext_pn, rmse_ext_tch, nf_ext)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_extended_residuals.csv"),
                    rows_raw, sel_extended_idx, res_ext)

    ctx_ho = build_context(rows_raw, sel_holdout_idx)
    res_ho, nf_ho = forward_apply(sel_holdout_idx, ctx_ho, st)
    rmse_ho_adn = rmse_pp(rows_raw, sel_holdout_idx, res_ho, :ADN)
    rmse_ho_pn  = rmse_pp(rows_raw, sel_holdout_idx, res_ho, :PN)
    rmse_ho_tch = rmse_pp(rows_raw, sel_holdout_idx, res_ho, :TCH)
    @printf("Holdout  RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp,  nfail=%d\n",
            rmse_ho_adn, rmse_ho_pn, rmse_ho_tch, nf_ho)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_holdout_residuals.csv"),
                    rows_raw, sel_holdout_idx, res_ho)

    # ---- Decision gates ----
    println("\n--- Decision gates ---")
    gate(name, val, thr) = begin
        ok = val < thr
        @printf("  [%s] %-42s  %.2f pp  (threshold < %.0f pp)\n",
                ok ? "PASS" : "FAIL", name, val, thr)
        return ok
    end

    g1 = gate("Core FE_ADN RMSE",      rmse_core_adn, 8.0)
    g2 = gate("Core FE_PN RMSE",       rmse_core_pn,  5.0)
    g3 = gate("Core FE_TCH RMSE",      rmse_core_tch, 4.0)
    g4 = gate("Extended FE_ADN RMSE",  rmse_ext_adn, 12.0)
    g5 = gate("Holdout FE_ADN RMSE",   rmse_ho_adn,  15.0)

    println("\n  All gates passed: ", all((g1, g2, g3, g4, g5)))
    println("=" ^ 72)
end

main()
