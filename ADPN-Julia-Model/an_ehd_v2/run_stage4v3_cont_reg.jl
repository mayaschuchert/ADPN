# -----------------------------------------------------------------------------
# run_stage4v3_cont_reg.jl — v3 joint LM with the new continuation solver,
# REGULARIZED toward the preserved (pre-kink-fix) θ.
#
# Hypothesis being tested: post-fix re-fits without regularization absorb
# φ_l-induced loss-surface bias (Core ADN gets better at the cost of Holdout
# generalization). Anchoring with L2 toward the preserved θ should hold the
# fit close to that point and only allow data-supported corrections.
#
# Output: output/stage4v3_cont_reg/
# -----------------------------------------------------------------------------
using Printf
using Dates
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.GalvContinuation

FitKinetics.SOLVER_BACKEND[] = :cont
@info "Solver backend set to" FitKinetics.SOLVER_BACKEND[]

# Preserved (pre-fix) θ — from output/preserved_pre_kink_fix/stage4_seq_fitted_theta.txt
const THETA_PRIOR = Float64[
    log10(6.505483e-3),    # log10 j0_ADN
    log10(1.212153e-3),    # log10 j0_PN
    log10(1.743811e-3),    # log10 j0_TCH
    0.536841,              # alpha_c_ADN
    0.521603,              # alpha_c_PN
    0.525030,              # alpha_c_TCH
    1.004517,              # n_ADN
    0.500000,              # n_PN  (will be clamped up to v3's LB=0.5; OK)
    1.000000,              # n_TCH (at v3's LB=1.0; OK)
]

# Per-parameter weights. Picked so each prior-residual term contributes
# roughly the same loss as ~1 pp data residual when θ moves 1 unit. Stronger
# weight on params that drifted to bounds without regularization (αc, n).
# This is a regularization knob — tune if needed.
const PRIOR_WEIGHT = Float64[
    1.0,    # log10 j0_ADN  — 1 pp² per decade
    1.0,    # log10 j0_PN
    1.0,    # log10 j0_TCH
    100.0,  # αc_ADN — strong: 1 pp² per 0.1 αc shift
    100.0,  # αc_PN
    100.0,  # αc_TCH
    25.0,   # n_ADN  — moderate
    25.0,   # n_PN
    25.0,   # n_TCH
]

# ---------- Paths ----------
const DATA_FILE    = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const OUT_DIR      = joinpath(@__DIR__, "output", "stage4v3_cont_reg")
const OUT_DATA_DIR = joinpath(OUT_DIR, "data")
const OUT_LOG_DIR  = joinpath(OUT_DIR, "logs")
isdir(OUT_DATA_DIR) || mkpath(OUT_DATA_DIR)
isdir(OUT_LOG_DIR)  || mkpath(OUT_LOG_DIR)

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

function write_residuals(path::String, rows::Vector{BloomquistRow},
                         sel::Vector{Int}, F::Vector{Float64})
    open(path, "w") do io
        println(io, "table,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN," *
                    "FE_ADN_obs,FE_ADN_model,FE_ADN_resid_pp," *
                    "FE_PN_obs,FE_PN_model,FE_PN_resid_pp," *
                    "FE_TCH_obs,FE_TCH_model,FE_TCH_resid_pp")
        for (n, idx) in pairs(sel)
            r = rows[idx]
            d_adn = F[3n - 2]; d_pn = F[3n - 1]; d_tch = F[3n]
            @printf(io, "%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                    r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                    r.FE_ADN_pct, r.FE_ADN_pct + d_adn, d_adn,
                    r.FE_PN_pct,  r.FE_PN_pct  + d_pn,  d_pn,
                    r.FE_TCH_pct, r.FE_TCH_pct + d_tch, d_tch)
        end
    end
end

function rmse(F::Vector{Float64}, ch::Symbol)
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

function _pin_flag(val::Float64, lb::Float64, ub::Float64; tol::Float64 = 1e-3)
    span = ub - lb
    (val - lb) ≤ tol * span && return "@LB"
    (ub - val) ≤ tol * span && return "@UB"
    return "  "
end

function main()
    println("="^72)
    println(" Stage 4v3-cont-REG — joint LM, continuation solver, L2 prior")
    println(" $(now())")
    println("="^72)

    rows_raw = load_bloomquist(DATA_FILE)
    sel_core_idx     = select_core(rows_raw)
    sel_extended_idx = select_extended(rows_raw)
    sel_holdout_idx  = select_holdout(rows_raw)
    @printf("Core / Extended / Holdout: %d / %d / %d rows\n",
            length(sel_core_idx), length(sel_extended_idx), length(sel_holdout_idx))

    println("\n--- Prior θ (preserved pre-kink-fix) ---")
    @printf("  j0_ADN=%.3e  j0_PN=%.3e  j0_TCH=%.3e\n",
            10^THETA_PRIOR[1], 10^THETA_PRIOR[2], 10^THETA_PRIOR[3])
    @printf("  αc:  ADN=%.3f  PN=%.3f  TCH=%.3f\n",
            THETA_PRIOR[4], THETA_PRIOR[5], THETA_PRIOR[6])
    @printf("  n:   ADN=%.3f  PN=%.3f  TCH=%.3f\n",
            THETA_PRIOR[7], THETA_PRIOR[8], THETA_PRIOR[9])
    @printf("  weights: %s\n", PRIOR_WEIGHT)

    println("\n--- Stage 4a: Regularized LM on Core ---")
    ctx_core = build_context(rows_raw, sel_core_idx)

    # Use the preserved θ as the starting point
    theta_init = clamp.(copy(THETA_PRIOR), THETA_LB, THETA_UB)
    @printf("theta_init = %s\n", theta_init)

    result = lm_fit(theta_init, ctx_core;
                    theta_prior  = THETA_PRIOR,
                    prior_weight = PRIOR_WEIGHT,
                    verbose      = true)
    println("\n[LM] done.  ", result.note)

    j0, ac, n_orders = theta_to_physical(result.theta)
    println("\nFitted parameters:")
    @printf("  j0_ADN = %.3e A/m²  %s    αc = %.3f  %s    n = %.3f  %s\n",
            j0[1], _pin_flag(result.theta[1], THETA_LB[1], THETA_UB[1]),
            ac[1], _pin_flag(result.theta[4], THETA_LB[4], THETA_UB[4]),
            n_orders[1], _pin_flag(result.theta[7], THETA_LB[7], THETA_UB[7]))
    @printf("  j0_PN  = %.3e A/m²  %s    αc = %.3f  %s    n = %.3f  %s\n",
            j0[2], _pin_flag(result.theta[2], THETA_LB[2], THETA_UB[2]),
            ac[2], _pin_flag(result.theta[5], THETA_LB[5], THETA_UB[5]),
            n_orders[2], _pin_flag(result.theta[8], THETA_LB[8], THETA_UB[8]))
    @printf("  j0_TCH = %.3e A/m²  %s    αc = %.3f  %s    n = %.3f  %s\n",
            j0[4], _pin_flag(result.theta[3], THETA_LB[3], THETA_UB[3]),
            ac[4], _pin_flag(result.theta[6], THETA_LB[6], THETA_UB[6]),
            n_orders[3], _pin_flag(result.theta[9], THETA_LB[9], THETA_UB[9]))

    # Compute Core RMSE — use plain residuals! (no regularization tail)
    F_core = zeros(3 * length(sel_core_idx))
    residuals!(F_core, result.theta, ctx_core)
    rmse_core_adn = rmse(F_core, :ADN)
    rmse_core_pn  = rmse(F_core, :PN)
    rmse_core_tch = rmse(F_core, :TCH)
    @printf("Core RMSE (data only) — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp\n",
            rmse_core_adn, rmse_core_pn, rmse_core_tch)

    # Save residuals + fitted θ
    write_residuals(joinpath(OUT_DATA_DIR, "stage4a_core_residuals.csv"),
                    rows_raw, sel_core_idx, F_core)

    open(joinpath(OUT_DATA_DIR, "stage4a_fitted_theta.txt"), "w") do io
        @printf(io, "# Stage 4v3-cont-REG — joint LM with L2 prior on preserved θ\n")
        @printf(io, "# date: %s\n", now())
        @printf(io, "# converged: %s   loss: %.6e\n", result.converged, result.loss)
        @printf(io, "# Core RMSE FE_ADN: %.4f pp,  FE_PN: %.4f pp,  FE_TCH: %.4f pp\n",
                rmse_core_adn, rmse_core_pn, rmse_core_tch)
        @printf(io, "# Prior weights: %s\n", PRIOR_WEIGHT)
        @printf(io, "j0_ADN   = %.6e\n", j0[1])
        @printf(io, "j0_PN    = %.6e\n", j0[2])
        @printf(io, "j0_HER   = %.6e  # frozen\n", J0_3_FROZEN)
        @printf(io, "j0_TCH   = %.6e\n", j0[4])
        @printf(io, "alpha_c_ADN = %.6f\n", ac[1])
        @printf(io, "alpha_c_PN  = %.6f\n", ac[2])
        @printf(io, "alpha_c_HER = %.6f  # frozen\n", ALPHA_C3_FROZEN)
        @printf(io, "alpha_c_TCH = %.6f\n", ac[4])
        @printf(io, "n_ADN = %.6f\n", n_orders[1])
        @printf(io, "n_PN  = %.6f\n", n_orders[2])
        @printf(io, "n_TCH = %.6f\n", n_orders[3])
    end

    # Stage 4b: Forward apply
    println("\n--- Stage 4b: forward apply to Extended + Holdout ---")
    ctx_ext = build_context(rows_raw, sel_extended_idx;
                            warm_init = copy(ctx_core.warm_by_key))
    F_ext = zeros(3 * length(ctx_ext.sel))
    residuals!(F_ext, result.theta, ctx_ext)
    rmse_ext_adn = rmse(F_ext, :ADN)
    rmse_ext_pn  = rmse(F_ext, :PN)
    rmse_ext_tch = rmse(F_ext, :TCH)
    @printf("Extended RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp\n",
            rmse_ext_adn, rmse_ext_pn, rmse_ext_tch)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_extended_residuals.csv"),
                    rows_raw, sel_extended_idx, F_ext)

    ctx_ho = build_context(rows_raw, sel_holdout_idx)
    F_ho = zeros(3 * length(ctx_ho.sel))
    residuals!(F_ho, result.theta, ctx_ho)
    rmse_ho_adn = rmse(F_ho, :ADN)
    rmse_ho_pn  = rmse(F_ho, :PN)
    rmse_ho_tch = rmse(F_ho, :TCH)
    @printf("Holdout  RMSE — FE_ADN: %.2f pp,  FE_PN: %.2f pp,  FE_TCH: %.2f pp\n",
            rmse_ho_adn, rmse_ho_pn, rmse_ho_tch)
    write_residuals(joinpath(OUT_DATA_DIR, "stage4b_holdout_residuals.csv"),
                    rows_raw, sel_holdout_idx, F_ho)

    println("\n--- Decision gates ---")
    gate(name, val, thr) = begin
        ok = val < thr
        @printf("  [%s] %-30s  %.2f pp  (< %.0f pp)\n",
                ok ? "PASS" : "FAIL", name, val, thr)
        return ok
    end
    g1 = gate("Core FE_ADN RMSE",     rmse_core_adn, 8.0)
    g2 = gate("Core FE_PN RMSE",      rmse_core_pn,  5.0)
    g3 = gate("Core FE_TCH RMSE",     rmse_core_tch, 4.0)
    g4 = gate("Extended FE_ADN RMSE", rmse_ext_adn, 12.0)
    g5 = gate("Holdout FE_ADN RMSE",  rmse_ho_adn,  15.0)
    println("\n  All gates passed: ", all((g1, g2, g3, g4, g5)))
    println("="^72)
end

main()
