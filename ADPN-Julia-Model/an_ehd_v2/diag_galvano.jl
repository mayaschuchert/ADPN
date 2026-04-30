
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.FitKinetics
using .ADPN_EHD.Chemistry
using .ADPN_EHD.FixedJ
using .ADPN_EHD.Mesh
using .ADPN_EHD.Assembly
using .ADPN_EHD.Kinetics
using .ADPN_EHD.Solver
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

function run_diag(gap_mm, Q, phi_AN, j_mA, theta_in, c_eq)
    j0, ac, n_orders = theta_to_physical(theta_in)
    Kinetics.set_kinetic_override!(j0, ac, (2.0, 1.0, 3.0))

    gap_m   = gap_mm * 1e-3
    Qt_m3s  = ADPN_EHD.Hydrodynamics.ml_min_to_m3_s(Q)
    delta_m = ADPN_EHD.Hydrodynamics.delta_leveque(gap_m, Qt_m3s)
    mesh    = ADPN_EHD.Mesh.make_mesh(100, delta_m; stretch = 10.0)
    u0      = make_initial_guess(100, c_eq, phi_AN)
    j_target = j_mA * 10.0

    @printf("\n=== gap=%.2f mm Q=%.0f phi=%.2f j_target=%.0f A/m² delta=%.0f µm ===\n",
            gap_mm, Q, phi_AN, j_target, delta_m * 1e6)

    full_res!(F, x, V) = full_residual!(F, x, mesh, phi_AN, V, 1.0, 1.0, c_eq)

    V_walk = -0.8; u_walk = copy(u0)
    V_prev_w = V_walk; u_prev_w = copy(u0)
    ds = 0.05; j_prev = 0.0
    stall_V = NaN; stall_u = nothing; stall_j = NaN

    for step in 1:300
        V_next = V_walk - ds
        if V_next < -2.5; break; end

        u_trial = copy(u_walk)
        res = newton_solve!(u_trial, (F, x) -> full_res!(F, x, V_next);
                            max_iter = 60, jacobian_mode = :fd, verbose = false)
        if !res.converged
            ds *= 0.3
            if ds < 5e-4
                @printf("STALL: V_curr=%.5f j_curr=%.2f A/m² (target=%.1f)\n",
                        V_walk, j_prev, j_target)
                stall_V = V_walk; stall_u = copy(u_walk); stall_j = j_prev
                break
            end
            continue
        end

        j1, j2, j3, j4 = faradaic_currents_from_state(u_trial, V_next, 1.0)
        j_now = j1 + j2 + j3 + j4

        V_prev_w = V_walk; u_prev_w = copy(u_walk)
        V_walk = V_next; u_walk = copy(u_trial); j_prev = j_now

        if res.iter <= 4;    ds = min(ds * 1.4, 0.20)
        elseif res.iter <= 10; ds = min(ds * 1.1, 0.20)
        elseif res.iter > 20;  ds = max(ds * 0.7, 5e-4)
        end

        if step % 20 == 0 || step <= 3
            @printf("  step=%3d V=%.4f j=%.2f A/m² iter=%d ds=%.4f\n",
                    step, V_walk, j_now, res.iter, ds)
        end

        if j_now >= j_target
            @printf("  Crossed at V=%.4f, j=%.2f\n", V_walk, j_now)
            stall_V = NaN
            break
        end
    end

    if !isnan(stall_V)
        F_check = zeros(length(stall_u))
        full_res!(F_check, stall_u, stall_V)
        @printf("normF at stall: max=%.3e, j_stall=%.2f, g=%.2f A/m²\n",
                maximum(abs.(F_check)), stall_j, stall_j - j_target)

        @printf("\n--- Galvanostatic Newton (verbose) ---\n")
        u_gal = copy(stall_u)
        build_res = (Vv) -> (F, x) -> full_residual!(F, x, mesh, phi_AN, Vv, 1.0, 1.0, c_eq)
        j_fn = (x, Vv) -> begin
            j1, j2, j3, j4 = faradaic_currents_from_state(x, Vv, 1.0)
            j1 + j2 + j3 + j4
        end
        res_g = newton_galvanostatic!(u_gal, stall_V, j_target, build_res, j_fn;
                                      tol_phys = 1e-5, tol_rel_j = 1e-3,
                                      max_iter = 30, max_backtrack = 15,
                                      verbose = true)
        @printf("Galvanostatic result: converged=%s iter=%d V=%.4f\n",
                res_g.converged, res_g.iter, res_g.V)
    end

    Kinetics.KIN_OVERRIDE[] = Kinetics.KIN_OVERRIDE[]  # restore not needed; just doc
end

rows_raw = load_bloomquist(DATA_FILE)
theta = clamp.(Float64[
    log10(6.505483e-3), log10(1.212153e-3), log10(1.743811e-3),
    0.536841, 0.521603, 0.525030,
    1.004517, 0.500000, 1.000000,
], THETA_LB, THETA_UB)
c_eq = Chemistry.solve_phosphate_equilibrium()

run_diag(0.50, 2.0, 0.16, 85.0, theta, c_eq)
run_diag(0.50, 2.0, 0.28, 177.0, theta, c_eq)
run_diag(1.00, 2.0, 0.28, 177.0, theta, c_eq)
