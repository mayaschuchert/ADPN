# Smoke test: solve one Bloomquist row with both solvers and compare.
using Printf
include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.GalvContinuation
using .ADPN_EHD.FixedJ

# Pick a benign row: gap=0.5mm, Q=2 mL/min, j=80 mA/cm², φ=0.29 (S5 first row)
gap_m  = 0.5e-3
Q_m3s  = ADPN_EHD.ml_min_to_m3_s(2.0)
delta  = ADPN_EHD.delta_leveque(gap_m, Q_m3s)
phi    = 0.29
j_tgt  = 80.0 * 10.0   # 80 mA/cm² → 800 A/m²

println("Test row: gap=0.5mm, Q=2 mL/min, φ=0.29, j_target=80 mA/cm²")
@printf("δ_lev = %.2f μm\n", delta * 1e6)

mesh = ADPN_EHD.make_mesh(100, delta; stretch = 10.0)
c_eq = ADPN_EHD.solve_phosphate_equilibrium()

# Use v6.x converged kinetics from stage4_seq fit (the preserved baseline)
j0  = (6.505483e-3, 1.212153e-3, 2.666e-5, 1.743811e-3)
ac  = (0.536841, 0.521603, 0.390, 0.525030)
n3  = (1.004517, 0.500000, 1.000000)

println("\n--- Old solver (V-walk + bisection) ---")
u_warm = ADPN_EHD.make_initial_guess(100, c_eq, phi)
t1 = @elapsed res_old = solve_at_j(j_tgt, phi, delta, mesh, u_warm, c_eq;
                                    j0 = j0, alpha_c = ac, n_orders = n3,
                                    V_hi = -0.3)
@printf("  converged = %s\n", res_old.converged)
@printf("  V = %.4f V,  j = %.2f A/m²,  FE_ADN = %.2f%%, FE_PN = %.2f%%, FE_HER = %.2f%%, FE_TCH = %.2f%%\n",
        res_old.V_cathode, res_old.j_total, res_old.FE_ADN_pct,
        res_old.FE_PN_pct, res_old.FE_HER_pct, res_old.FE_TCH_pct)
@printf("  wall = %.2f s, note = %s\n", t1, res_old.note)

println("\n--- New solver (galv-continuation) ---")
t2 = @elapsed res_new = solve_at_j_cont(j_tgt, phi, mesh, c_eq;
                                        j0 = j0, alpha_c = ac, n_orders = n3,
                                        verbose = false)
@printf("  converged = %s\n", res_new.converged)
@printf("  V = %.4f V,  j = %.2f A/m²,  FE_ADN = %.2f%%, FE_PN = %.2f%%, FE_HER = %.2f%%, FE_TCH = %.2f%%\n",
        res_new.V_cathode, res_new.j_total, res_new.FE_ADN_pct,
        res_new.FE_PN_pct, res_new.FE_HER_pct, res_new.FE_TCH_pct)
@printf("  wall = %.2f s, n_newton = %d, note = %s\n", t2, res_new.n_newton, res_new.note)

println("\n--- Comparison ---")
@printf("  ΔV       = %.4f V\n",   res_new.V_cathode - res_old.V_cathode)
@printf("  ΔFE_ADN  = %+.2f pp\n", res_new.FE_ADN_pct - res_old.FE_ADN_pct)
@printf("  ΔFE_PN   = %+.2f pp\n", res_new.FE_PN_pct  - res_old.FE_PN_pct)
@printf("  ΔFE_HER  = %+.2f pp\n", res_new.FE_HER_pct - res_old.FE_HER_pct)
@printf("  ΔFE_TCH  = %+.2f pp\n", res_new.FE_TCH_pct - res_old.FE_TCH_pct)

# ---- Test fixed_j clamping ----
println("\n--- Test FIXED_J clamping (HER fixed to 100 A/m²) ---")
j_her_clamp = 100.0   # A/m²
fixed = (NaN, NaN, j_her_clamp, NaN)
t3 = @elapsed res_clamped = solve_at_j_cont(j_tgt, phi, mesh, c_eq;
                                              j0 = j0, alpha_c = ac, n_orders = n3,
                                              fixed_j = fixed,
                                              verbose = false)
@printf("  converged = %s\n", res_clamped.converged)
@printf("  V = %.4f V,  j_total = %.2f A/m²\n", res_clamped.V_cathode, res_clamped.j_total)
@printf("  j_HER_actual = %.2f A/m²  (target = %.2f, error = %.4f)\n",
        res_clamped.j3, j_her_clamp, res_clamped.j3 - j_her_clamp)
@printf("  wall = %.2f s, n_newton = %d, note = %s\n", t3, res_clamped.n_newton, res_clamped.note)
