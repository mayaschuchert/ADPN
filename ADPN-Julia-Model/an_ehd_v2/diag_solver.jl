# Diagnostic: test solve_at_j on a few failing and succeeding rows with verbose=true
# to find the exact failure note.

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.Chemistry
using Printf

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")

using DelimitedFiles

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

rows_raw = load_bloomquist(DATA_FILE)
sel_core = select_core(rows_raw)
ctx = build_context(rows_raw, sel_core)

# Use stage4_seq theta as test point
theta = clamp.(Float64[
    log10(6.505483e-3), log10(1.212153e-3), log10(1.743811e-3),
    0.536841, 0.521603, 0.525030,
    1.004517, 0.500000, 1.000000,
], THETA_LB, THETA_UB)

j0, ac, n_orders = theta_to_physical(theta)

using .ADPN_EHD.FixedJ
using .ADPN_EHD.Mesh

# Test a sample of rows: 3 failing, 2 succeeding
test_cases = [
    ("S5 fail",  0.50, 2.0,  0.16, 850.0),   # Table_S5, Q=2, phi=0.16, j=85
    ("S5 fail2", 0.50, 2.0,  0.29, 1170.0),  # Table_S5, Q=2, phi=0.29, j=117
    ("S5 succ",  0.50, 2.0,  0.29, 800.0),   # Table_S5, Q=2, phi=0.29, j=80 → succeeds
    ("S6 succ",  0.50, 6.0,  0.29, 800.0),   # Table_S6, Q=6, phi=0.29, j=80 → succeeds
    ("S8 fail",  1.00, 2.0,  0.29, 800.0),   # Table_S8, Q=2, phi=0.29, j=80 → fails
]

c_eq = Chemistry.solve_phosphate_equilibrium()

for (label, gap_mm, Q, phi_AN, j_target) in test_cases
    gap_m   = gap_mm * 1e-3
    Qt_m3s  = ADPN_EHD.Hydrodynamics.ml_min_to_m3_s(Q)
    delta_m = ADPN_EHD.Hydrodynamics.delta_leveque(gap_m, Qt_m3s)
    delta_um = round(Int, delta_m * 1e6)
    mesh     = ADPN_EHD.Mesh.make_mesh(100, delta_m; stretch = 10.0)
    u0       = make_initial_guess(100, c_eq, phi_AN)

    @printf("\n=== %s  gap=%.1f mm  Q=%.0f  phi=%.2f  delta=%d µm  j=%.0f A/m² ===\n",
            label, gap_mm, Q, phi_AN, delta_um, j_target)

    for nmesh in (100, 200, 400)
        mesh_n = ADPN_EHD.Mesh.make_mesh(nmesh, delta_m; stretch = 10.0)
        u_n    = make_initial_guess(nmesh, c_eq, phi_AN)
        res = solve_at_j(j_target, phi_AN, delta_m, mesh_n, u_n, c_eq;
                         j0 = j0, alpha_c = ac, n_orders = n_orders,
                         verbose = false)
        @printf("  N=%d: converged=%s  note=\"%s\"\n",
                nmesh, res.converged, res.note)
        if res.converged
            @printf("    FE_ADN=%.1f%%  V=%.3f\n", res.FE_ADN_pct, res.V_cathode)
            break
        end
    end
end
