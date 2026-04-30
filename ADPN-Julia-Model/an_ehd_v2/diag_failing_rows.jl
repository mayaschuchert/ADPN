
# Identify which core rows still fail at the stage4_seq theta
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.Chemistry
using .ADPN_EHD.FixedJ
using .ADPN_EHD.Hydrodynamics
using Printf, DelimitedFiles

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")

function load_bloomquist(path)
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

rows_raw = load_bloomquist(DATA_FILE)
sel_core = select_core(rows_raw)
ctx = build_context(rows_raw, sel_core)

theta = clamp.(Float64[
    log10(6.505483e-3), log10(1.212153e-3), log10(1.743811e-3),
    0.536841, 0.521603, 0.525030,
    1.004517, 0.500000, 1.000000,
], THETA_LB, THETA_UB)

j0, ac, n_orders = theta_to_physical(theta)
c_eq = Chemistry.solve_phosphate_equilibrium()

@printf("\nChecking all %d core rows at theta_init:\n\n", length(sel_core))
@printf("%-10s %6s %6s %6s %6s %8s %8s %12s\n",
        "table", "gap_mm", "Q", "phi_AN", "j_mA", "delta_µm", "j_lim_mA", "status")
println("-"^80)

using .ADPN_EHD.Params
using .ADPN_EHD.Kinetics

global nfail = 0
for (i, idx) in enumerate(sel_core)
    row = rows_raw[idx]
    gap_m   = row.gap_mm * 1e-3
    Qt_m3s  = ml_min_to_m3_s(row.Q_total_mL_min)
    delta_m = delta_leveque(gap_m, Qt_m3s)
    mesh    = ADPN_EHD.Mesh.make_mesh(100, delta_m; stretch = 10.0)
    u0      = make_initial_guess(100, c_eq, row.phi_AN)
    j_target = row.j_mA_cm2 * 10.0   # A/m²

    # Rough transport limit for AN (n_e=2, D_AN=2.30e-9 m²/s)
    # c_AN_bulk in mol/m³ from phi_AN and C_AN_SAT
    c_AN_bulk = row.phi_AN * ADPN_EHD.Params.C_AN_SAT   # mol/m³ (rough)
    j_lim = 2 * 96485.0 * 2.30e-9 * c_AN_bulk / delta_m / 10.0  # mA/cm²

    res = solve_at_j(j_target, row.phi_AN, delta_m, mesh, u0, c_eq;
                     j0 = j0, alpha_c = ac, n_orders = n_orders,
                     verbose = false)

    status = res.converged ? "OK($(res.note[1:min(4,end)]))" : "FAIL"
    if !res.converged
        global nfail += 1
        @printf("%-10s %6.2f %6.1f %6.3f %6.1f %8.1f %8.1f %12s  <-- FAIL note=\"%s\"\n",
                row.table, row.gap_mm, row.Q_total_mL_min, row.phi_AN,
                row.j_mA_cm2, delta_m*1e6, j_lim, status, res.note)
    end
end
@printf("\n%d/%d core rows fail at theta_init\n", nfail, length(sel_core))
