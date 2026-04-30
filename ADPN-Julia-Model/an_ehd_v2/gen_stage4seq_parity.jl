# gen_stage4seq_parity.jl — run stage4_seq fitted kinetics on all rows and emit
# a CSV that the parity-plot Python script reads.

using Printf
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const THETA_FILE = joinpath(@__DIR__, "output", "stage4_seq", "data",
                            "stage4_seq_fitted_theta.txt")
const OUT_CSV   = joinpath(@__DIR__, "output", "stage4_seq", "data",
                           "stage4seq_parity.csv")

# ---------- load fitted theta from text file ----------
function read_theta_file(path)
    vals = Dict{String,Float64}()
    for line in readlines(path)
        m = match(r"^(\w+)\s*=\s*([\d.e+\-]+)", line)
        m === nothing && continue
        vals[m[1]] = parse(Float64, m[2])
    end
    return vals
end

# ---------- load rows (same as run_stage4_seq) ----------
function load_bloomquist(path)
    raw, hdr = readdlm(path, ','; header = true)
    hdr = vec(hdr)
    col(n) = findfirst(==(string(n)), hdr)
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

println("Loading data...")
rows_raw   = load_bloomquist(DATA_FILE)
sel_core   = select_core(rows_raw)
sel_ext    = select_extended(rows_raw)
sel_hold   = select_holdout(rows_raw)

# set label per row index (core overrides extended since sel_core ⊂ sel_ext)
set_label = fill("other", length(rows_raw))
for i in sel_ext;    set_label[i] = "extended";  end
for i in sel_hold;   set_label[i] = "holdout";   end
for i in sel_core;   set_label[i] = "core";      end   # must be last

@printf("core=%d  extended=%d  holdout=%d  total=%d\n",
        length(sel_core), length(sel_ext), length(sel_hold), length(rows_raw))

# ---------- read fitted kinetics ----------
tv = read_theta_file(THETA_FILE)
j0_fit = (tv["j0_ADN"],      tv["j0_PN"],       tv["j0_HER"],      tv["j0_TCH"])
ac_fit = (tv["alpha_c_ADN"], tv["alpha_c_PN"],   tv["alpha_c_HER"], tv["alpha_c_TCH"])
n_fit  = (tv["n_ADN"],       tv["n_PN"],         tv["n_TCH"])

@printf("Fitted: j0=(%.3e, %.3e, %.3e, %.3e)\n", j0_fit...)
@printf("        ac=(%.3f, %.3f, %.3f, %.3f)\n",  ac_fit...)
@printf("        n=(%.3f, %.3f, %.3f)\n",          n_fit...)
flush(stdout)

# ---------- build context for all unique conditions ----------
all_sel = union(sel_core, sel_ext, sel_hold) |> collect |> sort
ctx = build_context(rows_raw, all_sel)
sol_cache = Dict{Any, Tuple{Float64, Vector{Float64}}}()

# ---------- evaluate every row ----------
header = ["row_idx","set","gap_mm","Q_total_mL_min","j_mA_cm2","phi_AN",
          "We_aq","We_org",
          "FE_ADN_obs","FE_PN_obs","FE_TCH_obs",
          "FE_ADN_pred","FE_PN_pred","FE_TCH_pred",
          "V_cathode_pred","converged"]

results = Any[]
nfail = 0

for idx in all_sel
    r    = ctx.rows[idx]          # use ctx.rows so delta_lev_m is populated
    key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
    δkey = round(r.delta_lev_m * 1e6; digits = 1)
    mesh = ctx.mesh_by_delta[δkey]
    u_w  = get(ctx.warm_by_key, key,
                make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN))

    V_h, u_h = get(sol_cache, key, (NaN, nothing))

    result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                        mesh, u_w, ctx.c_eq;
                        j0 = j0_fit, alpha_c = ac_fit, n_orders = n_fit,
                        V_hint = V_h, u_hint = u_h,
                        ds_min = 1e-5, newton_max_iter = 120)

    if !result.converged && occursin("already ≥ j_target", result.note)
        result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                            mesh, u_w, ctx.c_eq;
                            j0 = j0_fit, alpha_c = ac_fit, n_orders = n_fit,
                            ds_min = 1e-5, newton_max_iter = 120)
    end

    if result.converged
        sol_cache[key] = (result.V_cathode, result.state)
        ctx.warm_by_key[key] = result.state
        push!(results, [idx, set_label[idx],
                        r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.We_aq, r.We_org,
                        r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct,
                        result.FE_ADN_pct, result.FE_PN_pct, result.FE_TCH_pct,
                        result.V_cathode, 1])
        @printf("  row %3d (%s)  V=%.3f  FE_ADN: obs=%.1f pred=%.1f\n",
                idx, set_label[idx], result.V_cathode,
                r.FE_ADN_pct, result.FE_ADN_pct)
    else
        global nfail += 1
        push!(results, [idx, set_label[idx],
                        r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN,
                        r.We_aq, r.We_org,
                        r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct,
                        NaN, NaN, NaN, NaN, 0])
        @printf("  row %3d (%s)  FAIL: %s\n", idx, set_label[idx], result.note)
    end
    flush(stdout)
end

@printf("\nDone: %d rows  nfail=%d\n", length(results), nfail)

open(OUT_CSV, "w") do io
    println(io, join(header, ","))
    for row in results
        println(io, join(row, ","))
    end
end
println("Saved: ", OUT_CSV)
