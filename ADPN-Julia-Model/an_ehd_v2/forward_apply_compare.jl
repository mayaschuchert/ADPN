# forward_apply_compare.jl — run two θ candidates through the new continuation
# solver on all 162 Bloomquist rows, writing parity-format CSVs that the
# existing plot_weber_maps_model.py can read.
#
# Outputs:
#   output/forward_v3/data/parity.csv          — v3 fit θ (joint LM, V-walk-trained)
#   output/forward_stage4seq_LB/data/parity.csv — stage4_seq overwritten θ
#                                                (n_PN=0.500@LB, n_TCH=1.000@LB)

using Printf
using Dates
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.GalvContinuation

FitKinetics.SOLVER_BACKEND[] = :cont
@info "Solver backend" FitKinetics.SOLVER_BACKEND[]

# v3 θ — from output/stage4v3/data/stage4a_fitted_theta.txt
const THETA_V3 = (
    j0_ADN  = 8.454694e-3,  ac_ADN = 0.566594, n_ADN = 1.000000,
    j0_PN   = 1.616525e-3,  ac_PN  = 0.560814, n_PN  = 1.636977,
    j0_HER  = 2.666e-5,     ac_HER = 0.390,
    j0_TCH  = 2.426481e-3,  ac_TCH = 0.563606, n_TCH = 1.549037,
)

# stage4_seq LB-pinned θ — from output/stage4_seq/data/stage4_seq_fitted_theta.txt
# (the version that overwrote at 10:28, NOT the recovered Run-2 θ)
const THETA_SEQ_LB = (
    j0_ADN  = 6.505483e-3,  ac_ADN = 0.536841, n_ADN = 1.004517,
    j0_PN   = 1.212153e-3,  ac_PN  = 0.521603, n_PN  = 0.500000,
    j0_HER  = 2.666e-5,     ac_HER = 0.390,
    j0_TCH  = 1.743811e-3,  ac_TCH = 0.525030, n_TCH = 1.000000,
)

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")

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

function evaluate_set(ctx, sel, t)
    j0_t = (t.j0_ADN, t.j0_PN, t.j0_HER, t.j0_TCH)
    ac_t = (t.ac_ADN, t.ac_PN, t.ac_HER, t.ac_TCH)
    n_t  = (t.n_ADN, t.n_PN, t.n_TCH)
    n = length(sel)
    Vcat       = fill(NaN, n)
    FEpre      = fill((NaN, NaN, NaN, NaN), n)
    converged  = falses(n)
    nfail = 0
    for (i, idx) in pairs(sel)
        r = ctx.rows[idx]
        key = round(r.delta_lev_m * 1e6; digits = 1)
        mesh = ctx.mesh_by_delta[key]
        result = solve_at_j_cont(r.j_target_A_m2, r.phi_AN, mesh, ctx.c_eq;
                                  j0 = j0_t, alpha_c = ac_t, n_orders = n_t,
                                  verbose = false)
        if result.converged
            Vcat[i]  = result.V_cathode
            FEpre[i] = (result.FE_ADN_pct, result.FE_PN_pct,
                         result.FE_HER_pct, result.FE_TCH_pct)
            converged[i] = true
        else
            nfail += 1
        end
    end
    return Vcat, FEpre, converged, nfail
end

function dump_parity(out_csv::String, rows, sel_core, sel_ext, sel_ho,
                      ctx_core, ctx_ext, ctx_ho, t)
    Vc, Fc, Cc, nfc = evaluate_set(ctx_core, sel_core, t)
    Ve, Fe, Ce, nfe = evaluate_set(ctx_ext,  sel_ext,  t)
    Vh, Fh, Ch, nfh = evaluate_set(ctx_ho,   sel_ho,   t)
    @printf("    Core nfail=%d/%d, Ext nfail=%d/%d, Holdout nfail=%d/%d\n",
            nfc, length(sel_core), nfe, length(sel_ext), nfh, length(sel_ho))

    open(out_csv, "w") do io
        println(io, "row_idx,set,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN,We_aq,We_org," *
                    "FE_ADN_obs,FE_PN_obs,FE_TCH_obs," *
                    "FE_ADN_pred,FE_PN_pred,FE_HER_pred,FE_TCH_pred," *
                    "V_cathode_pred,converged")
        function dump_set(set_name, sel, V, F, C)
            for (i, idx) in pairs(sel)
                r = rows[idx]
                fA, fP, fH, fT = F[i]
                @printf(io, "%d,%s,%.2f,%.1f,%.0f,%.4f,%.4f,%.4e,%.2f,%.2f,%.2f,",
                        idx, set_name, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2,
                        r.phi_AN, r.We_aq, r.We_org,
                        r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct)
                if C[i]
                    @printf(io, "%.2f,%.2f,%.2f,%.2f,%.4f,1\n",
                            fA, fP, fH, fT, V[i])
                else
                    @printf(io, "NA,NA,NA,NA,NA,0\n")
                end
            end
        end
        dump_set("core",     sel_core, Vc, Fc, Cc)
        dump_set("extended", sel_ext,  Ve, Fe, Ce)
        dump_set("holdout",  sel_ho,   Vh, Fh, Ch)
    end
end

function main()
    rows_raw = load_bloomquist(DATA_FILE)
    sel_core = select_core(rows_raw)
    sel_ext  = select_extended(rows_raw)
    sel_ho   = select_holdout(rows_raw)
    @printf("Loaded %d rows.  Core=%d Ext=%d Ho=%d\n",
            length(rows_raw), length(sel_core), length(sel_ext), length(sel_ho))

    ctx_core = build_context(rows_raw, sel_core)
    ctx_ext  = build_context(rows_raw, sel_ext)
    ctx_ho   = build_context(rows_raw, sel_ho)

    cases = [
        ("forward_v3",            "v3 θ (joint LM, post-fix)",       THETA_V3),
        ("forward_stage4seq_LB",  "stage4_seq LB-pinned θ (post-fix)", THETA_SEQ_LB),
    ]

    for (subdir, label, t) in cases
        outdir = joinpath(@__DIR__, "output", subdir, "data")
        isdir(outdir) || mkpath(outdir)
        @printf("\n=== %s ===\n", label)
        @printf("  Writing → %s/parity.csv\n", outdir)
        dump_parity(joinpath(outdir, "parity.csv"),
                    rows_raw, sel_core, sel_ext, sel_ho,
                    ctx_core, ctx_ext, ctx_ho, t)
    end
end

main()
