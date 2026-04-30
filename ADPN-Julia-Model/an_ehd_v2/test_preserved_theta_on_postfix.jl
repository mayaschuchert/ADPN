# Forward-apply the preserved (pre-kink-fix) fitted θ to the post-kink-fix
# physics. No fitting — just one solve per row at the saved parameters.
#
# Question this answers: did the preserved kinetic params absorb the Na+/φ_l
# bug (so that running them on the corrected physics gives garbage), or do
# they still produce competitive RMSEs under the corrected physics?

using Printf
using Dates
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.GalvContinuation

# ---- Switch to the new continuation solver ----
FitKinetics.SOLVER_BACKEND[] = :cont

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const OUT_DIR   = joinpath(@__DIR__, "output", "preserved_theta_check")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---- Preserved θ from stage4_seq_fitted_theta.txt (pre-kink-fix run) ----
const J0_ADN_PRE = 6.505483e-3
const J0_PN_PRE  = 1.212153e-3
const J0_HER_PRE = 2.666e-5         # frozen
const J0_TCH_PRE = 1.743811e-3
const AC_ADN_PRE = 0.536841
const AC_PN_PRE  = 0.521603
const AC_HER_PRE = 0.390000         # frozen
const AC_TCH_PRE = 0.525030
const N_ADN_PRE  = 1.004517
const N_PN_PRE   = 0.500000         # @LB (script's LB=0.5)
const N_TCH_PRE  = 1.000000         # @LB (script's LB=1.0)

# ---- CSV loader (same as run_stage4v3.jl) ----
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
println(" Forward-apply preserved θ (pre-kink-fix) on post-kink-fix physics")
println(" $(now())")
println("="^72)

rows_raw = load_bloomquist(DATA_FILE)
sel_core_idx     = select_core(rows_raw)
sel_extended_idx = select_extended(rows_raw)
sel_holdout_idx  = select_holdout(rows_raw)
@printf("Core / Extended / Holdout: %d / %d / %d rows\n",
        length(sel_core_idx), length(sel_extended_idx), length(sel_holdout_idx))

# Build a v3-style theta vector to feed into the residuals! function.
# v3 layout: (log10 j0_ADN, log10 j0_PN, log10 j0_TCH, αc_ADN, αc_PN, αc_TCH,
#             n_ADN, n_PN, n_TCH).  HER frozen via constants in FitKinetics.
theta = Float64[
    log10(J0_ADN_PRE),  log10(J0_PN_PRE),  log10(J0_TCH_PRE),
    AC_ADN_PRE,         AC_PN_PRE,         AC_TCH_PRE,
    N_ADN_PRE,          N_PN_PRE,          N_TCH_PRE
]

# v3's THETA_LB has n_ADN ≥ 1.0 and n_PN ≥ 0.5; the preserved values are at or
# above those bounds, so clamping to v3 bounds is a no-op for them.
theta_clamped = clamp.(theta, THETA_LB, THETA_UB)
all_in_bounds = theta_clamped == theta
@printf("\nPreserved θ within v3 bounds: %s\n", all_in_bounds)
if !all_in_bounds
    println("  Differences (preserved → clamped):")
    for k in 1:9
        if theta[k] != theta_clamped[k]
            @printf("    [%d] %s → %s\n", k, theta[k], theta_clamped[k])
        end
    end
end

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

function evaluate(rows, sel, name::String)
    println("\n--- $name (n=$(length(sel))) ---")
    ctx = build_context(rows, sel)
    F = zeros(3 * length(sel))
    nfail = residuals!(F, theta, ctx)   # uses cont solver via SOLVER_BACKEND[]
    rA = rmse_pp(F, :ADN)
    rP = rmse_pp(F, :PN)
    rT = rmse_pp(F, :TCH)
    @printf("  RMSE FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp,  nfail=%d/%d\n",
            rA, rP, rT, nfail, length(sel))
    return (rA, rP, rT, nfail, F)
end

cA, cP, cT, cN, cF = evaluate(rows_raw, sel_core_idx,     "Core")
eA, eP, eT, eN, _  = evaluate(rows_raw, sel_extended_idx, "Extended")
hA, hP, hT, hN, _  = evaluate(rows_raw, sel_holdout_idx,  "Holdout")

println("\n--- Comparison vs v3 (Core gates) ---")
println("                              v3       preserved-on-postfix")
@printf("  Core FE_ADN RMSE:        %.2f pp     %.2f pp\n", 8.44, cA)
@printf("  Core FE_PN  RMSE:        %.2f pp     %.2f pp\n", 5.32, cP)
@printf("  Core FE_TCH RMSE:        %.2f pp     %.2f pp\n", 4.55, cT)
@printf("  Holdout ADN RMSE:       %.2f pp    %.2f pp\n",  30.68, hA)

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

# Save Core residuals so we can compare row-by-row to v3 if needed
open(joinpath(OUT_DIR, "core_residuals_preserved_theta.csv"), "w") do io
    println(io, "table,gap_mm,Q,j_mA_cm2,phi_AN,FE_ADN_obs,FE_ADN_resid_pp," *
                "FE_PN_obs,FE_PN_resid_pp,FE_TCH_obs,FE_TCH_resid_pp")
    for (n, idx) in pairs(sel_core_idx)
        r = rows_raw[idx]
        @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                r.FE_ADN_pct, cF[3n-2],
                r.FE_PN_pct,  cF[3n-1],
                r.FE_TCH_pct, cF[3n  ])
    end
end
println("\nWrote per-row Core residuals to $(OUT_DIR)/core_residuals_preserved_theta.csv")
println("="^72)
