# -----------------------------------------------------------------------------
# run_stage4.jl — DEPRECATED. v6 6-parameter fit. Outputs in output/stage4/
# are the canonical v6 baseline and must not be overwritten.
#
# Use run_stage4v2.jl for the v6.x 8-parameter fit (adds n_1, n_2 — AN
# reaction orders for ADPN and PN). v2 writes to output/stage4v2/.
#
# This script will refuse to run; remove the error() below if you really
# need to re-execute the original 6-param fit (you'll need to also revert
# fit_kinetics.jl to N_THETA=6 and the corresponding bounds/THETA0 entries).
# -----------------------------------------------------------------------------
error("run_stage4.jl is frozen as the v6 baseline. Use run_stage4v2.jl " *
      "for the v6.x fit with n_1, n_2 freed; outputs under output/stage4v2/.")

# -----------------------------------------------------------------------------
# (original v6 6-param pipeline kept below for reference)
# -----------------------------------------------------------------------------
using Printf
using Dates
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics

# ---------- Paths ----------
const DATA_FILE     = joinpath(@__DIR__, "Experimental_data", "bloomquist_data.csv")
const OUT_DIR       = joinpath(@__DIR__, "output", "stage4")
const OUT_DATA_DIR  = joinpath(OUT_DIR, "data")
const OUT_LOG_DIR   = joinpath(OUT_DIR, "logs")
isdir(OUT_DATA_DIR) || mkpath(OUT_DATA_DIR)
isdir(OUT_LOG_DIR)  || mkpath(OUT_LOG_DIR)

# ---------- Minimal CSV loader (Base only — no DataFrames dependency) ----------
function load_bloomquist(path::String)
    raw, hdr = readdlm(path, ','; header = true)
    hdr = vec(hdr)
    col(name) = findfirst(==(string(name)), hdr)
    @assert col("table") !== nothing "missing column: table"
    @assert col("FE_ADN_pct") !== nothing "missing column: FE_ADN_pct"

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
            # Derived fields filled by build_context:
            NaN, NaN, NaN, NaN
        ))
    end
    return rows
end

# ---------- Residual CSV writer ----------
function write_residuals(path::String, rows::Vector{BloomquistRow},
                         sel::Vector{Int}, F::Vector{Float64},
                         theta::Vector{Float64})
    open(path, "w") do io
        println(io, "table,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN," *
                    "FE_ADN_obs,FE_ADN_model,FE_ADN_resid_pp," *
                    "FE_PN_obs,FE_PN_model,FE_PN_resid_pp")
        for (n, idx) in pairs(sel)
            r = rows[idx]
            d_adn = F[2n - 1]
            d_pn  = F[2n    ]
            mod_adn = r.FE_ADN_pct + d_adn
            mod_pn  = r.FE_PN_pct  + d_pn
            @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                    r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                    r.FE_ADN_pct, mod_adn, d_adn,
                    r.FE_PN_pct,  mod_pn,  d_pn)
        end
    end
end

function rmse(F::Vector{Float64}, channel::Symbol)
    n = length(F) ÷ 2
    s = 0.0
    if channel === :ADN
        for i in 1:n; s += F[2i - 1]^2; end
    elseif channel === :PN
        for i in 1:n; s += F[2i    ]^2; end
    else
        error("channel must be :ADN or :PN")
    end
    return sqrt(s / n)
end

# ---------- Main ----------
function main()
    println("=" ^ 72)
    println(" Stage 4 — Bloomquist kinetics-only fit (v6 §20)")
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

    # ---------- Stage 4a — fit on Core ----------
    println("\n--- Stage 4a: LM fit on Core ---")
    ctx_core = build_context(rows_raw, sel_core_idx)
    @printf("Pre-cached %d unique-δ meshes\n", length(ctx_core.mesh_by_delta))

    result = lm_fit(THETA0, ctx_core; verbose = true)
    println("\n[LM] done.  ", result.note)

    j0, ac = theta_to_physical(result.theta)
    @printf("  j0_1 = %.3e A/m²    α_c1 = %.3f\n", j0[1], ac[1])
    @printf("  j0_2 = %.3e A/m²    α_c2 = %.3f\n", j0[2], ac[2])
    @printf("  j0_3 = %.3e A/m²    α_c3 = %.3f\n", j0[3], ac[3])

    # Final residuals at fitted θ on Core (for CSV + RMSE)
    F_core = zeros(2 * length(ctx_core.sel))
    residuals!(F_core, result.theta, ctx_core)
    rmse_core_adn = rmse(F_core, :ADN)
    rmse_core_pn  = rmse(F_core, :PN)
    @printf("Core RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp\n",
            rmse_core_adn, rmse_core_pn)

    write_residuals(joinpath(OUT_DATA_DIR, "stage4a_core_residuals.csv"),
                    rows_raw, sel_core_idx, F_core, result.theta)

    # Save fitted parameters
    open(joinpath(OUT_DATA_DIR, "stage4a_fitted_theta.txt"), "w") do io
        @printf(io, "# v6 Stage 4a — fitted kinetic parameters (transport frozen)\n")
        @printf(io, "# date: %s\n", now())
        @printf(io, "# converged: %s\n", result.converged)
        @printf(io, "# loss: %.6e (sum of squared residuals, pp²)\n", result.loss)
        @printf(io, "# Core rows: %d\n", length(sel_core_idx))
        @printf(io, "# RMSE FE_ADN: %.4f pp\n", rmse_core_adn)
        @printf(io, "# RMSE FE_PN:  %.4f pp\n", rmse_core_pn)
        @printf(io, "j0_1 = %.6e\n", j0[1])
        @printf(io, "j0_2 = %.6e\n", j0[2])
        @printf(io, "j0_3 = %.6e\n", j0[3])
        @printf(io, "alpha_c1 = %.6f\n", ac[1])
        @printf(io, "alpha_c2 = %.6f\n", ac[2])
        @printf(io, "alpha_c3 = %.6f\n", ac[3])
    end

    # ---------- Stage 4b — forward apply, no re-fit ----------
    println("\n--- Stage 4b: forward apply Stage 4a θ to Extended + Holdout ---")

    # Re-use the warm cache from Core where keys overlap; build new ctxs
    ctx_ext = build_context(rows_raw, sel_extended_idx;
                             warm_init = copy(ctx_core.warm_by_key))
    F_ext = zeros(2 * length(ctx_ext.sel))
    residuals!(F_ext, result.theta, ctx_ext)
    rmse_ext_adn = rmse(F_ext, :ADN)
    rmse_ext_pn  = rmse(F_ext, :PN)
    @printf("Extended RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp\n",
            rmse_ext_adn, rmse_ext_pn)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_extended_residuals.csv"),
                    rows_raw, sel_extended_idx, F_ext, result.theta)

    ctx_ho = build_context(rows_raw, sel_holdout_idx)
    F_ho = zeros(2 * length(ctx_ho.sel))
    residuals!(F_ho, result.theta, ctx_ho)
    rmse_ho_adn = rmse(F_ho, :ADN)
    rmse_ho_pn  = rmse(F_ho, :PN)
    @printf("Holdout  RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp\n",
            rmse_ho_adn, rmse_ho_pn)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_holdout_residuals.csv"),
                    rows_raw, sel_holdout_idx, F_ho, result.theta)

    # ---------- Decision-gate summary (§20.4) ----------
    println("\n--- Decision gates (§20.4) ---")
    gate(name, val, thr, dir = :lt) = begin
        ok = dir === :lt ? val < thr : val > thr
        marker = ok ? "PASS" : "FAIL"
        @printf("  [%s] %-40s  %.2f pp  (threshold: %s %.0f pp)\n",
                marker, name, val, dir === :lt ? "<" : ">", thr)
        return ok
    end

    g1 = gate("Core FE_ADN RMSE",       rmse_core_adn, 8.0)
    g2 = gate("Core FE_PN RMSE",        rmse_core_pn,  5.0)
    g3 = gate("Extended FE_ADN RMSE",   rmse_ext_adn, 12.0)
    g4 = gate("Holdout FE_ADN RMSE",    rmse_ho_adn,  15.0)

    println("\n  All gates passed: ", all((g1, g2, g3, g4)))
    println("=" ^ 72)
end

main()
