# -----------------------------------------------------------------------------
# analyze_stage4.jl — post-fit diagnostic.
#
# Reads the fitted kinetics from output/stage4/data/stage4a_fitted_theta.txt,
# re-runs the model at each Bloomquist row, and emits an enriched residual
# CSV with the back-derived experimental V_cell:
#
#     V_cell_obs = PR_ADN [kg/cm²/h]  /  (EP_ADN [kg/kWh] · j [A/cm²])
#
# alongside the model's V_cell_pred = V_CE + |V_cathode| + j · R_series.
# Useful for evaluating whether the frozen v6 defaults V_CE = 1.7 V and
# R_contact = 1e-4 Ω·m² are reasonable, and as a v7 input for fitting them.
#
# Run AFTER run_stage4.jl has produced stage4a_fitted_theta.txt.
# -----------------------------------------------------------------------------
using Printf
using Dates
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD

const DATA_FILE      = joinpath(@__DIR__, "Experimental_data", "bloomquist_data.csv")
const STAGE4_DIR     = joinpath(@__DIR__, "output", "stage4", "data")
const FITTED_THETA   = joinpath(STAGE4_DIR, "stage4a_fitted_theta.txt")
const OUT_CSV        = joinpath(STAGE4_DIR, "stage4_diagnostic.csv")

# ---------- Reuse the loader from run_stage4.jl by inlining ----------
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

# Parse `key = value` lines from stage4a_fitted_theta.txt
function load_fitted_theta(path::String)
    vals = Dict{String,Float64}()
    for line in eachline(path)
        startswith(strip(line), "#") && continue
        m = match(r"^\s*(\w+)\s*=\s*([\d\-\+\.eE]+)", line)
        m === nothing && continue
        vals[m.captures[1]] = parse(Float64, m.captures[2])
    end
    j0 = (vals["j0_1"], vals["j0_2"], vals["j0_3"])
    ac = (vals["alpha_c1"], vals["alpha_c2"], vals["alpha_c3"])
    return j0, ac
end

# Back-derive V_cell [V] from EP/PR/j (NaN where degenerate).
#
# Dimensional check:
#   EP [kg/kWh] = mass / (P [kW] · t [h])
#   PR [kg/cm²/h] = mass / (A [cm²] · t [h])
#   ⇒ P [kW] = A · PR / EP                                       [cm² · kg/cm²/h ÷ kg/kWh]
#   P [W] = j [A/cm²] · A [cm²] · V_cell [V]
#   ⇒ V_cell [V] = P [W] / (j · A)
#                = 1000 · A · PR / EP / (j · A)
#                = 1000 · PR / (EP · j)
@inline function v_cell_obs(r::BloomquistRow)
    j_A_cm2 = r.j_mA_cm2 * 1e-3
    if r.PR_ADN_kg_cm2_h <= 1e-10 || r.EP_ADN_kg_kWh <= 1e-6 || j_A_cm2 <= 0
        return NaN
    end
    return 1000.0 * r.PR_ADN_kg_cm2_h / (r.EP_ADN_kg_kWh * j_A_cm2)
end

# Subset tag (Core / Extended-only / Holdout / Excluded) for one row index.
function subset_tag(r::BloomquistRow)
    in_core = (r.gap_mm == 0.5 || r.gap_mm == 1.0) &&
              r.j_mA_cm2 ≤ 190.0 &&
              r.phi_AN ≥ 0.04
    in_ext  = (r.gap_mm == 0.5 || r.gap_mm == 1.0) && r.phi_AN ≥ 0.04
    in_ho   = r.gap_mm == 0.25 && r.phi_AN ≥ 0.04
    if in_core
        return "Core"
    elseif in_ext
        return "Extended-only"
    elseif in_ho
        return "Holdout"
    else
        return "Excluded"
    end
end

function main()
    println("=" ^ 72)
    println(" Post-fit V_cell diagnostic — ", now())
    println("=" ^ 72)

    j0, ac = load_fitted_theta(FITTED_THETA)
    @printf("Loaded fitted kinetics from %s\n", basename(FITTED_THETA))
    @printf("  j0_1 = %.3e A/m²    α_c1 = %.3f\n", j0[1], ac[1])
    @printf("  j0_2 = %.3e A/m²    α_c2 = %.3f\n", j0[2], ac[2])
    @printf("  j0_3 = %.3e A/m²    α_c3 = %.3f\n", j0[3], ac[3])

    rows = load_bloomquist(DATA_FILE)
    @printf("Loaded %d rows from %s\n", length(rows), basename(DATA_FILE))

    # Build a context over ALL rows so we get a full warm-start cache and one
    # mesh per unique δ.
    sel_all = collect(1:length(rows))
    ctx = build_context(rows, sel_all)
    @printf("Pre-cached %d unique-δ meshes\n\n", length(ctx.mesh_by_delta))

    c_eq = ctx.c_eq
    open(OUT_CSV, "w") do io
        println(io, "table,subset,gap_mm,Q_total_mL_min,j_mA_cm2,phi_AN," *
                    "FE_ADN_obs,FE_ADN_model,FE_ADN_resid_pp," *
                    "FE_PN_obs,FE_PN_model,FE_PN_resid_pp," *
                    "V_cathode_SHE_V,V_cell_obs_V,V_cell_pred_V,V_cell_resid_V," *
                    "delta_lev_um,kappa_eff_S_per_m,R_series_Ohm_m2," *
                    "converged,note")

        sums = Dict("Core"=>Float64[], "Extended-only"=>Float64[],
                    "Holdout"=>Float64[], "Excluded"=>Float64[])
        nrows_processed = 0

        for (idx, r) in pairs(ctx.rows)
            tag = subset_tag(r)
            mesh_k = round(r.delta_lev_m * 1e6; digits = 1)
            mesh = get(ctx.mesh_by_delta, mesh_k, nothing)
            if mesh === nothing
                # Build on demand for any δ we missed
                mesh = make_mesh(ctx.N_mesh, r.delta_lev_m;
                                  stretch = ctx.stretch)
                ctx.mesh_by_delta[mesh_k] = mesh
            end
            wkey = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
            u_warm = if haskey(ctx.warm_by_key, wkey)
                         ctx.warm_by_key[wkey]
                     else
                         make_initial_guess(ctx.N_mesh, c_eq, r.phi_AN)
                     end

            res = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                             mesh, u_warm, c_eq;
                             j0 = j0, alpha_c = ac, verbose = false)

            if res.converged
                ctx.warm_by_key[wkey] = res.state
            end

            v_obs = v_cell_obs(r)
            v_pred = res.converged ?
                CellVoltage.V_cell_predicted(res.V_cathode, r.j_target_A_m2,
                                             r.gap_m, r.delta_lev_m,
                                             r.phi_AN, c_eq) : NaN
            v_resid = (isfinite(v_obs) && isfinite(v_pred)) ?
                      (v_pred - v_obs) : NaN
            κeff = CellVoltage.kappa_eff(c_eq, r.phi_AN)

            d_adn = res.converged ? res.FE_ADN_pct - r.FE_ADN_pct : NaN
            d_pn  = res.converged ? res.FE_PN_pct  - r.FE_PN_pct  : NaN

            @printf(io,
                "%s,%s,%.2f,%.0f,%.0f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.2f,%.3f,%.4e,%s,%s\n",
                r.table, tag, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                r.FE_ADN_pct,
                res.converged ? res.FE_ADN_pct : NaN, d_adn,
                r.FE_PN_pct,
                res.converged ? res.FE_PN_pct : NaN, d_pn,
                res.converged ? res.V_cathode : NaN,
                v_obs, v_pred, v_resid,
                r.delta_lev_m * 1e6, κeff, r.R_series_Ohm_m2,
                res.converged, replace(res.note, "," => ";"))

            if isfinite(v_resid)
                push!(sums[tag], v_resid)
            end
            nrows_processed += 1
            if nrows_processed % 20 == 0
                @printf("  %d/%d rows processed...\n",
                        nrows_processed, length(ctx.rows))
            end
        end

        println()
        for (tag, name) in [("Core","Core"),("Extended-only","Extended-only"),
                            ("Holdout","Holdout"),("Excluded","Excluded")]
            arr = sums[tag]
            n   = length(arr)
            if n == 0
                @printf("  %-15s :  (no rows with valid V_cell_obs)\n", name)
                continue
            end
            med  = sort(arr)[(n + 1) ÷ 2]
            mae  = sum(abs, arr) / n
            bias = sum(arr) / n
            @printf("  %-15s :  n=%3d   median(resid)=%+.3f V   MAE=%.3f V   bias=%+.3f V\n",
                    name, n, med, mae, bias)
        end
    end
    @printf("\nWrote %s\n", OUT_CSV)
    println("=" ^ 72)
end

main()
