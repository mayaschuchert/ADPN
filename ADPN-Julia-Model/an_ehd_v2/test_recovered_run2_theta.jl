# Forward-apply the RECOVERED Run-2 stage4_seq θ (with interior n values, all
# Core gates passing, ~707 pp² loss) on the post-kink-fix physics.
#
# These params are stored at:
#   output/preserved_pre_kink_fix/stage4_seq/data/stage4_seq_fitted_theta_RUN2_recovered.txt
#
# Question: are they usable under the corrected physics (Na+ in φ_l)?

using Printf
using Dates
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.GalvContinuation

FitKinetics.SOLVER_BACKEND[] = :cont

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const OUT_DIR   = joinpath(@__DIR__, "output", "recovered_run2_check")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# Recovered Run-2 params — interior n's, all Core gates passed pre-fix
const J0_ADN  = 2.703e-3
const J0_PN   = 1.152e-3
const J0_HER  = 2.666e-5      # frozen
const J0_TCH  = 1.701e-3
const AC_ADN  = 0.567
const AC_PN   = 0.507
const AC_HER  = 0.390         # frozen
const AC_TCH  = 0.519
const N_ADN   = 0.629
const N_PN    = 0.966
const N_TCH   = 1.826

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

println("="^72)
println(" Forward-apply RECOVERED Run-2 θ on post-kink-fix physics")
println(" $(now())")
println("="^72)
@printf("  j0_ADN=%.3e  αc_ADN=%.3f  n_ADN=%.3f\n", J0_ADN, AC_ADN, N_ADN)
@printf("  j0_PN =%.3e  αc_PN =%.3f  n_PN =%.3f\n", J0_PN,  AC_PN,  N_PN)
@printf("  j0_HER=%.3e  αc_HER=%.3f  (frozen)\n",  J0_HER, AC_HER)
@printf("  j0_TCH=%.3e  αc_TCH=%.3f  n_TCH=%.3f\n", J0_TCH, AC_TCH, N_TCH)

rows_raw = load_bloomquist(DATA_FILE)
sel_core = select_core(rows_raw)
sel_ext  = select_extended(rows_raw)
sel_ho   = select_holdout(rows_raw)
@printf("\nCore / Extended / Holdout: %d / %d / %d rows\n",
        length(sel_core), length(sel_ext), length(sel_ho))

# v3-style theta vector — recovered Run 2 had n_ADN=0.629 below v3's LB=1.0.
# Don't clamp to v3 bounds; bypass `residuals!` and call solver directly so we
# can keep n_ADN at its original 0.629 value (post-fix model's physics doesn't
# care about the LM bound — that's only used during fitting).
function rmse_pp(F::Vector{Float64}, ch::Symbol)
    n = length(F) ÷ 3
    s = 0.0
    if ch === :ADN
        for i in 1:n; s += F[3i-2]^2; end
    elseif ch === :PN
        for i in 1:n; s += F[3i-1]^2; end
    elseif ch === :TCH
        for i in 1:n; s += F[3i  ]^2; end
    end
    return sqrt(s / n)
end

# Use solve_at_j_cont directly to bypass any LM-bound clamping
function evaluate_set(rows, sel, name::String)
    println("\n--- $name (n=$(length(sel))) ---")
    ctx = build_context(rows, sel)
    F = zeros(3 * length(sel))
    nfail = 0
    j0_t = (J0_ADN, J0_PN, J0_HER, J0_TCH)
    ac_t = (AC_ADN, AC_PN, AC_HER, AC_TCH)
    n_t  = (N_ADN, N_PN, N_TCH)
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
        else
            F[3i-2] = 1000.0
            F[3i-1] = 1000.0
            F[3i  ] = 1000.0
            nfail += 1
        end
    end
    rA = rmse_pp(F, :ADN); rP = rmse_pp(F, :PN); rT = rmse_pp(F, :TCH)
    @printf("  RMSE FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp,  nfail=%d/%d\n",
            rA, rP, rT, nfail, length(sel))
    return rA, rP, rT, nfail, F
end

cA, cP, cT, cN, cF = evaluate_set(rows_raw, sel_core, "Core")
eA, eP, eT, eN, _  = evaluate_set(rows_raw, sel_ext,  "Extended")
hA, hP, hT, hN, _  = evaluate_set(rows_raw, sel_ho,   "Holdout")

println("\n" * "─"^72)
println("  Comparison vs other candidates")
println("─"^72)
@printf("                              Core ADN  Core PN  Core TCH  Holdout ADN\n")
@printf("  Run-2 recovered θ:           %5.2f    %5.2f    %5.2f      %5.2f\n", cA, cP, cT, hA)
@printf("  Preserved (LB-pinned) θ:     13.86    4.88     2.83       25.99\n")
@printf("  v3 (post-fix re-fit):         8.44    5.32     4.55       30.68\n")
@printf("  v3-cont (post-fix re-fit):    9.07    5.42     3.41       29.97\n")
@printf("  stage4_seq_cont (re-fit):     9.77    4.32     3.83       26.65\n")

println("\n--- Decision gates ---")
gate(name, val, thr) = begin
    ok = val < thr
    @printf("  [%s] %-30s  %.2f pp  (< %.0f pp)\n",
            ok ? "PASS" : "FAIL", name, val, thr)
    return ok
end
g1 = gate("Core FE_ADN RMSE",     cA,  8.0)
g2 = gate("Core FE_PN RMSE",      cP,  5.0)
g3 = gate("Core FE_TCH RMSE",     cT,  4.0)
g4 = gate("Extended FE_ADN RMSE", eA, 12.0)
g5 = gate("Holdout FE_ADN RMSE",  hA, 15.0)
println("\n  All gates passed: ", all((g1, g2, g3, g4, g5)))

# Save Core residuals
open(joinpath(OUT_DIR, "core_residuals_recovered_run2.csv"), "w") do io
    println(io, "table,gap_mm,Q,j_mA_cm2,phi_AN,FE_ADN_obs,FE_ADN_resid_pp," *
                "FE_PN_obs,FE_PN_resid_pp,FE_TCH_obs,FE_TCH_resid_pp")
    for (n, idx) in pairs(sel_core)
        r = rows_raw[idx]
        @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                r.FE_ADN_pct, cF[3n-2],
                r.FE_PN_pct,  cF[3n-1],
                r.FE_TCH_pct, cF[3n  ])
    end
end
println("\nWrote per-row Core residuals to $(OUT_DIR)/core_residuals_recovered_run2.csv")
println("="^72)
