# Diagnose which Core rows fail under S1.11 clamping at the v4 converged θ.
# For each branch (HER → ADN → PN → TCH), call branch_residuals! once and
# log which row indices return non-converged solutions.

using Printf
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.FitKineticsSeq
using .ADPN_EHD.GalvContinuation
using .ADPN_EHD.Kinetics

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

println("Loading Bloomquist data + Core selection...")
rows_raw  = load_bloomquist(DATA_FILE)
sel_core  = select_core(rows_raw)
ctx_core  = build_context(rows_raw, sel_core)
@printf("Core: %d rows; %d unique-δ meshes\n",
        length(sel_core), length(ctx_core.mesh_by_delta))

# Build SeqFitState at the v4 converged θ
st = SeqFitState()
st.j0_ADN, st.ac_ADN, st.n_ADN = 1.473e-3, 0.563, 1.000
st.j0_PN,  st.ac_PN,  st.n_PN  = 1.279e-3, 0.498, 1.025
st.j0_HER, st.ac_HER           = 3.031e-5, 0.395
st.j0_TCH, st.ac_TCH, st.n_TCH = 9.470e-4, 0.566, 1.519

# theta_branch placeholder — we just want to evaluate at the committed values
# (re-passing them as the active branch parameters so they're applied)
function eval_branch(branch::Symbol)
    nres = length(sel_core)
    F = zeros(nres)
    theta = if branch === :HER
        [log10(st.j0_HER), st.ac_HER]
    elseif branch === :ADN
        [log10(st.j0_ADN), st.ac_ADN, st.n_ADN]
    elseif branch === :PN
        [log10(st.j0_PN),  st.ac_PN,  st.n_PN]
    elseif branch === :TCH
        [log10(st.j0_TCH), st.ac_TCH, st.n_TCH]
    end
    nfail = branch_residuals!(F, theta, st, ctx_core, branch)
    return F, nfail
end

println("\n" * "="^72)
for branch in (:HER, :ADN, :PN, :TCH)
    println("\n--- Branch: $branch ---")
    F, nfail = eval_branch(branch)
    @printf("  nfail = %d / %d rows\n", nfail, length(sel_core))
    @printf("  loss  = %.4e (log10(j+1) units²)\n", sum(F.^2))
    failed_idx = findall(x -> abs(x - 5.0) < 1e-6, F)   # 5.0 = the penalty
    if !isempty(failed_idx)
        println("  Failing rows (table, gap_mm, Q, j_mA_cm2, φ_AN, FE_ADN, FE_PN, FE_TCH, FE_HER_derived):")
        for i in failed_idx
            r = ctx_core.rows[sel_core[i]]
            fe_her = max(0.0, 100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct))
            @printf("    [row %3d] %s  gap=%.2fmm Q=%.0f j=%.0f φ=%.4f  FE_ADN=%.1f FE_PN=%.1f FE_TCH=%.1f FE_HER=%.1f\n",
                    sel_core[i], r.table, r.gap_mm, r.Q_total_mL_min,
                    r.j_mA_cm2, r.phi_AN,
                    r.FE_ADN_pct, r.FE_PN_pct, r.FE_TCH_pct, fe_her)
        end
    end
end

println("\n" * "="^72)
println("Summary: rows failing in ALL 4 branches")
all_failures = Set{Int}()
fail_per_branch = Dict{Symbol, Set{Int}}()
for branch in (:HER, :ADN, :PN, :TCH)
    F, _ = eval_branch(branch)
    failed_idx = findall(x -> abs(x - 5.0) < 1e-6, F)
    fail_set = Set(sel_core[i] for i in failed_idx)
    fail_per_branch[branch] = fail_set
    union!(all_failures, fail_set)
end
common_failures = intersect(values(fail_per_branch)...)
println("\nUnion of failing rows across branches: $(length(all_failures)) rows")
println("Intersection (failing in ALL 4 branches): $(length(common_failures)) rows")
for idx in sort(collect(common_failures))
    r = rows_raw[idx]
    fe_her = max(0.0, 100.0 - (r.FE_ADN_pct + r.FE_PN_pct + r.FE_TCH_pct))
    @printf("  [row %3d] %s  gap=%.2fmm Q=%.0f j=%.0f φ=%.4f  FE_HER_derived=%.1f%%\n",
            idx, r.table, r.gap_mm, r.Q_total_mL_min, r.j_mA_cm2, r.phi_AN, fe_her)
end
