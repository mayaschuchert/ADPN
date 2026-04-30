# publish_recovered_run2.jl — Final publication artifacts.
#
# Uses the RECOVERED Run-2 stage4_seq θ (interior n's, all Core gates passed
# pre-fix, recovered from output/stage4_seq/logs/stage4_seq.log line 347)
# as the published fit, but generates plot-ready data by FORWARD-APPLYING those
# parameters to the post-kink-fix physics (Na+ now in φ_l current conservation).
#
# Outputs (in output/publish_recovered_run2/):
#   data/fitted_theta.txt              — published kinetic params + provenance
#   data/parity.csv                    — row_idx,set,...,FE_*_obs,FE_*_pred,V_cathode,converged
#                                         (matches stage4seq_parity.csv format so existing
#                                         plot_stage4seq_*.py scripts work)
#   data/core_residuals.csv            — per-row Core residuals
#   data/extended_residuals.csv        — per-row Extended residuals
#   data/holdout_residuals.csv         — per-row Holdout residuals
#   data/summary.txt                   — RMSEs + gate decisions

using Printf
using Dates
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.GalvContinuation

FitKinetics.SOLVER_BACKEND[] = :cont

# ---- Recovered Run-2 θ (preserved at
#      output/preserved_pre_kink_fix/stage4_seq/data/stage4_seq_fitted_theta_RUN2_recovered.txt) ----
const J0_ADN  = 2.703e-3
const J0_PN   = 1.152e-3
const J0_HER  = 2.666e-5      # frozen (v6.x converged)
const J0_TCH  = 1.701e-3
const AC_ADN  = 0.567
const AC_PN   = 0.507
const AC_HER  = 0.390         # frozen (v6.x converged)
const AC_TCH  = 0.519
const N_ADN   = 0.629
const N_PN    = 0.966
const N_TCH   = 1.826

const DATA_FILE    = joinpath(@__DIR__, "Experimental_data", "bloomquist_data.csv")
const OUT_DIR      = joinpath(@__DIR__, "output")
const OUT_DATA_DIR = joinpath(OUT_DIR, "data")
isdir(OUT_DATA_DIR) || mkpath(OUT_DATA_DIR)

# ---------- CSV loader ----------
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

# ---------- Forward-apply ----------
function evaluate_set(ctx, sel)
    j0_t = (J0_ADN, J0_PN, J0_HER, J0_TCH)
    ac_t = (AC_ADN, AC_PN, AC_HER, AC_TCH)
    n_t  = (N_ADN, N_PN, N_TCH)
    n = length(sel)
    F = zeros(3 * n)
    Vcat   = fill(NaN, n)
    FEpre  = fill((NaN, NaN, NaN, NaN), n)
    converged = falses(n)
    nfail = 0

    for (i, idx) in pairs(sel)
        r = ctx.rows[idx]
        key = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[key]
        result = solve_at_j_cont(r.j_target_A_m2, r.phi_AN, mesh, ctx.c_eq;
                                  j0 = j0_t, alpha_c = ac_t, n_orders = n_t,
                                  verbose = false)
        if result.converged
            F[3i-2] = result.FE_ADN_pct - r.FE_ADN_pct
            F[3i-1] = result.FE_PN_pct  - r.FE_PN_pct
            F[3i  ] = result.FE_TCH_pct - r.FE_TCH_pct
            Vcat[i]  = result.V_cathode
            FEpre[i] = (result.FE_ADN_pct, result.FE_PN_pct,
                         result.FE_HER_pct, result.FE_TCH_pct)
            converged[i] = true
        else
            F[3i-2] = NaN; F[3i-1] = NaN; F[3i] = NaN
            nfail += 1
        end
    end
    return (F = F, Vcat = Vcat, FEpre = FEpre, converged = converged, nfail = nfail)
end

function rmse_pp(F::Vector{Float64}, ch::Symbol)
    n = length(F) ÷ 3
    s = 0.0; cnt = 0
    for i in 1:n
        v = ch === :ADN ? F[3i-2] : ch === :PN ? F[3i-1] : F[3i]
        if !isnan(v)
            s += v^2; cnt += 1
        end
    end
    return cnt == 0 ? NaN : sqrt(s / cnt)
end

# ---------- Writers ----------
function write_residuals(path::String, rows, sel, F::Vector{Float64})
    open(path, "w") do io
        println(io, "table,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN," *
                    "FE_ADN_obs,FE_ADN_model,FE_ADN_resid_pp," *
                    "FE_PN_obs,FE_PN_model,FE_PN_resid_pp," *
                    "FE_TCH_obs,FE_TCH_model,FE_TCH_resid_pp")
        for (i, idx) in pairs(sel)
            r = rows[idx]
            d_adn = F[3i-2]; d_pn = F[3i-1]; d_tch = F[3i]
            if isnan(d_adn)
                @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,NA,NA,%.2f,NA,NA,%.2f,NA,NA\n",
                        r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct)
            else
                @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                        r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.FE_ADN_pct, r.FE_ADN_pct + d_adn, d_adn,
                        r.FE_PN_pct,  r.FE_PN_pct  + d_pn,  d_pn,
                        r.FE_TCH_pct, r.FE_TCH_pct + d_tch, d_tch)
            end
        end
    end
end

# ---------- Main ----------
function main()
    println("="^72)
    println(" publish_recovered_run2.jl — Final fit artifacts")
    println(" $(now())")
    println("="^72)
    @printf("Published θ: j0_ADN=%.3e  j0_PN=%.3e  j0_TCH=%.3e\n", J0_ADN, J0_PN, J0_TCH)
    @printf("             αc:  ADN=%.3f  PN=%.3f  TCH=%.3f  (HER=%.3f frozen)\n",
            AC_ADN, AC_PN, AC_TCH, AC_HER)
    @printf("             n:   ADN=%.3f  PN=%.3f  TCH=%.3f\n", N_ADN, N_PN, N_TCH)

    rows_raw = load_bloomquist(DATA_FILE)
    sel_core = select_core(rows_raw)
    sel_ext  = select_extended(rows_raw)
    sel_ho   = select_holdout(rows_raw)
    @printf("Rows: Core=%d  Extended=%d  Holdout=%d\n",
            length(sel_core), length(sel_ext), length(sel_ho))

    println("\n--- Forward-apply on post-fix physics ---")
    ctx_core = build_context(rows_raw, sel_core)
    rc = evaluate_set(ctx_core, sel_core)
    @printf("Core      RMSE  ADN=%.2f  PN=%.2f  TCH=%.2f  nfail=%d/%d\n",
            rmse_pp(rc.F, :ADN), rmse_pp(rc.F, :PN), rmse_pp(rc.F, :TCH),
            rc.nfail, length(sel_core))

    ctx_ext = build_context(rows_raw, sel_ext;
                             warm_init = copy(ctx_core.warm_by_key))
    re = evaluate_set(ctx_ext, sel_ext)
    @printf("Extended  RMSE  ADN=%.2f  PN=%.2f  TCH=%.2f  nfail=%d/%d\n",
            rmse_pp(re.F, :ADN), rmse_pp(re.F, :PN), rmse_pp(re.F, :TCH),
            re.nfail, length(sel_ext))

    ctx_ho = build_context(rows_raw, sel_ho)
    rh = evaluate_set(ctx_ho, sel_ho)
    @printf("Holdout   RMSE  ADN=%.2f  PN=%.2f  TCH=%.2f  nfail=%d/%d\n",
            rmse_pp(rh.F, :ADN), rmse_pp(rh.F, :PN), rmse_pp(rh.F, :TCH),
            rh.nfail, length(sel_ho))

    # ---- fitted_theta.txt ----
    open(joinpath(OUT_DATA_DIR, "fitted_theta.txt"), "w") do io
        @printf(io, "# Final published fit — recovered Run-2 stage4_seq θ\n")
        @printf(io, "# date: %s\n", now())
        @printf(io, "# Provenance:\n")
        @printf(io, "#   - Original fit: stage4_seq.log line 347 (2026-04-29 ~02:30)\n")
        @printf(io, "#   - Recovered file: output/preserved_pre_kink_fix/stage4_seq/data/\n")
        @printf(io, "#                     stage4_seq_fitted_theta_RUN2_recovered.txt\n")
        @printf(io, "#   - Plot data: forward-applied to POST-kink-fix physics (Na+ in φ_l)\n")
        @printf(io, "#                via galvanostatic continuation solver\n")
        @printf(io, "#\n")
        @printf(io, "# Post-fix RMSEs (this run):\n")
        @printf(io, "#   Core      ADN=%.2f pp  PN=%.2f pp  TCH=%.2f pp\n",
                rmse_pp(rc.F, :ADN), rmse_pp(rc.F, :PN), rmse_pp(rc.F, :TCH))
        @printf(io, "#   Extended  ADN=%.2f pp  PN=%.2f pp  TCH=%.2f pp\n",
                rmse_pp(re.F, :ADN), rmse_pp(re.F, :PN), rmse_pp(re.F, :TCH))
        @printf(io, "#   Holdout   ADN=%.2f pp  PN=%.2f pp  TCH=%.2f pp\n",
                rmse_pp(rh.F, :ADN), rmse_pp(rh.F, :PN), rmse_pp(rh.F, :TCH))
        @printf(io, "j0_ADN      = %.6e\n", J0_ADN)
        @printf(io, "j0_PN       = %.6e\n", J0_PN)
        @printf(io, "j0_HER      = %.6e  # frozen (v6.x converged)\n", J0_HER)
        @printf(io, "j0_TCH      = %.6e\n", J0_TCH)
        @printf(io, "alpha_c_ADN = %.6f\n", AC_ADN)
        @printf(io, "alpha_c_PN  = %.6f\n", AC_PN)
        @printf(io, "alpha_c_HER = %.6f  # frozen (v6.x converged)\n", AC_HER)
        @printf(io, "alpha_c_TCH = %.6f\n", AC_TCH)
        @printf(io, "n_ADN       = %.6f\n", N_ADN)
        @printf(io, "n_PN        = %.6f\n", N_PN)
        @printf(io, "n_TCH       = %.6f\n", N_TCH)
    end

    # ---- per-row residuals (3 sets) ----
    write_residuals(joinpath(OUT_DATA_DIR, "core_residuals.csv"),
                    rows_raw, sel_core, rc.F)
    write_residuals(joinpath(OUT_DATA_DIR, "extended_residuals.csv"),
                    rows_raw, sel_ext, re.F)
    write_residuals(joinpath(OUT_DATA_DIR, "holdout_residuals.csv"),
                    rows_raw, sel_ho, rh.F)

    # ---- combined parity CSV in stage4seq_parity.csv format ----
    open(joinpath(OUT_DATA_DIR, "parity.csv"), "w") do io
        println(io, "row_idx,set,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN,We_aq,We_org," *
                    "FE_ADN_obs,FE_PN_obs,FE_TCH_obs," *
                    "FE_ADN_pred,FE_PN_pred,FE_HER_pred,FE_TCH_pred," *
                    "V_cathode_pred,converged")
        function dump_set(set_name, sel, res)
            for (i, idx) in pairs(sel)
                r = rows_raw[idx]
                fA, fP, fH, fT = res.FEpre[i]
                @printf(io, "%d,%s,%.2f,%.1f,%.0f,%.4f,%.4f,%.4e,%.2f,%.2f,%.2f,",
                        idx, set_name, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2,
                        r.phi_AN, r.We_aq, r.We_org,
                        r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct)
                if res.converged[i]
                    @printf(io, "%.2f,%.2f,%.2f,%.2f,%.4f,1\n",
                            fA, fP, fH, fT, res.Vcat[i])
                else
                    @printf(io, "NA,NA,NA,NA,NA,0\n")
                end
            end
        end
        dump_set("core",     sel_core, rc)
        dump_set("extended", sel_ext,  re)
        dump_set("holdout",  sel_ho,   rh)
    end

    # ---- summary ----
    open(joinpath(OUT_DATA_DIR, "summary.txt"), "w") do io
        @printf(io, "# Summary — recovered Run-2 θ on post-kink-fix physics\n")
        @printf(io, "# date: %s\n\n", now())
        @printf(io, "Set        ADN_RMSE    PN_RMSE    TCH_RMSE   nfail\n")
        @printf(io, "Core        %5.2f       %5.2f       %5.2f      %d/%d\n",
                rmse_pp(rc.F, :ADN), rmse_pp(rc.F, :PN), rmse_pp(rc.F, :TCH),
                rc.nfail, length(sel_core))
        @printf(io, "Extended    %5.2f       %5.2f       %5.2f      %d/%d\n",
                rmse_pp(re.F, :ADN), rmse_pp(re.F, :PN), rmse_pp(re.F, :TCH),
                re.nfail, length(sel_ext))
        @printf(io, "Holdout     %5.2f       %5.2f       %5.2f      %d/%d\n",
                rmse_pp(rh.F, :ADN), rmse_pp(rh.F, :PN), rmse_pp(rh.F, :TCH),
                rh.nfail, length(sel_ho))
        @printf(io, "\nDecision gates:\n")
        gate_str(name, val, thr) =
            @sprintf("  [%s] %-30s  %.2f pp  (< %.0f pp)\n",
                     val < thr ? "PASS" : "FAIL", name, val, thr)
        print(io, gate_str("Core FE_ADN RMSE",     rmse_pp(rc.F, :ADN),  8.0))
        print(io, gate_str("Core FE_PN RMSE",      rmse_pp(rc.F, :PN),   5.0))
        print(io, gate_str("Core FE_TCH RMSE",     rmse_pp(rc.F, :TCH),  4.0))
        print(io, gate_str("Extended FE_ADN RMSE", rmse_pp(re.F, :ADN), 12.0))
        print(io, gate_str("Holdout FE_ADN RMSE",  rmse_pp(rh.F, :ADN), 15.0))
    end

    println("\nWrote $(OUT_DATA_DIR):")
    for f in readdir(OUT_DATA_DIR)
        @printf("  %s\n", f)
    end
    println("="^72)
end

main()
