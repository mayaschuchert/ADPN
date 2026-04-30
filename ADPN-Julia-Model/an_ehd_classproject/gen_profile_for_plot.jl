# gen_profile_for_plot.jl  — generate a species profile CSV at a representative
# operating condition (phi_AN=0.05, gap=1mm, Q=50 mL/min, j≈100 A/m²)
# using Stage 3 cache warm start if available.

using Printf
using DelimitedFiles
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics

OUT = joinpath(@__DIR__, "output", "data", "profile_operating.csv")

rows_raw = begin
    raw, hdr = readdlm(joinpath(@__DIR__, "Experimental_data",
                                "bloomquist_data.csv"), ','; header=true)
    hdr = vec(hdr)
    col(n) = findfirst(==(string(n)), hdr)
    rows = BloomquistRow[]
    for i in 1:size(raw,1)
        push!(rows, BloomquistRow(
            string(raw[i,col("table")]),
            Float64(raw[i,col("gap_mm")]),
            Float64(raw[i,col("Q_total_mL_min")]),
            Float64(raw[i,col("j_mA_cm2")]),
            Float64(raw[i,col("phi_AN")]),
            Float64(raw[i,col("Q_aq_mL_min")]),
            Float64(raw[i,col("Q_org_mL_min")]),
            Float64(raw[i,col("We_aq")]),
            Float64(raw[i,col("We_org")]),
            Float64(raw[i,col("FE_ADN_pct")]),
            Float64(raw[i,col("FE_TCH_pct")]),
            Float64(raw[i,col("FE_PN_pct")]),
            Float64(raw[i,col("PR_ADN_kg_cm2_h")]),
            Float64(raw[i,col("EP_ADN_kg_kWh")]),
            NaN, NaN, NaN, NaN))
    end
    rows
end

# Pick a representative row from the core set (moderate current, mid phi_AN)
sel = select_core(rows_raw)
ctx = build_context(rows_raw, sel)

# Find a row with j_target near 100 A/m² and phi_AN near 0.05
target_idx = sel[argmin([abs(ctx.rows[i].j_target_A_m2 - 100.0) for i in sel])]
r = ctx.rows[target_idx]
@printf("Using row: gap=%.1f mm  Q=%.0f mL/min  phi_AN=%.3f  j_target=%.1f A/m²\n",
        r.gap_mm, r.Q_total_mL_min, r.phi_AN, r.j_target_A_m2)

δkey = round(r.delta_lev_m * 1e6; digits=1)
mesh = ctx.mesh_by_delta[δkey]
key  = (r.gap_mm, r.Q_total_mL_min, r.phi_AN)
u_w  = get(ctx.warm_by_key, key, make_initial_guess(ctx.N_mesh, ctx.c_eq, r.phi_AN))

j0_t, ac_t, n_t = FitKinetics.theta_to_physical(THETA0)
result = solve_at_j(r.j_target_A_m2, r.phi_AN, r.delta_lev_m,
                    mesh, u_w, ctx.c_eq;
                    j0=j0_t, alpha_c=ac_t, n_orders=n_t)

if !result.converged
    println("solve_at_j did not converge: ", result.note, " — using warm start")
    state = u_w
else
    @printf("Converged: V=%.4f V  j=%.1f A/m²  FE_ADN=%.1f%%  FE_PN=%.1f%%  FE_TCH=%.1f%%\n",
            result.V_cathode, result.j_total,
            result.FE_ADN_pct, result.FE_PN_pct, result.FE_TCH_pct)
    state = result.state
end

# Decode state into concentration matrix and phi_l
N = ctx.N_mesh
c, phi = decode_state(state, N)

# Accumulate x positions from mesh
x = cumsum([0.0; mesh.dx]) .- mesh.dx[1]/2
x = [mesh.dx[k]/2 + sum(mesh.dx[1:k-1]) for k in 1:N]

header = ["x_m","c_H","c_OH","c_H2PO4","c_HPO4","c_PO4",
          "c_AN","c_ADPN","c_PN","c_TCH","phi_l"]

data = hcat(x,
            c[1,:], c[2,:], c[3,:],
            c[4,:], c[5,:], c[6,:],
            c[7,:], c[8,:], c[9,:],
            phi)

open(OUT, "w") do io
    println(io, join(header, ","))
    writedlm(io, data, ',')
end
println("Saved: ", OUT)
