# Test: can log-j continuation proceed past V=-2.32 if we start from an
# earlier well-converged state (V=-2.0 V)?
using Printf

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.Params, .ADPN_EHD.Chemistry, .ADPN_EHD.Mesh
using .ADPN_EHD.Assembly, .ADPN_EHD.Solver, .ADPN_EHD.Kinetics

mesh = make_mesh(100, 50e-6; stretch = 10.0)
c_eq = solve_phosphate_equilibrium()
u = make_initial_guess(100, c_eq, 0.0)

# Bootstrap
for α in range(0.1, 1.0; length = 10)
    res_cb = (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, α, 0.0, c_eq)
    newton_solve!(u, res_cb)
end
for α in [min(1e-6 * 2.0^k, 1.0) for k in 0:20]
    res_cb = (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, 1.0, α, c_eq)
    newton_solve!(u, res_cb)
end

# V-continuation to V = -2.0
build_res = V -> (F, u) -> full_residual!(F, u, mesh, 0.0, V, 1.0, 1.0, c_eq)
hist = newton_continuation(u, -1.0, -2.0, build_res;
                           ds_init = 0.05, ds_min = 1e-4, max_iter = 40, verbose = false)
V_at, u_at = hist[end]
@printf("V-continuation to V=%.4f OK (%d points)\n", V_at, length(hist))

# Now log-j from V=-2.0 to V=-2.5
faradaic_current = (u_state, V_now) -> begin
    c_AN_s = exp(clamp(u_state[6], -50.0, 50.0))
    phi_s  = u_state[9]
    j1, j2, j3 = tafel_currents(c_AN_s, phi_s, V_now, 1.0)
    return j1 + j2 + j3
end

j0 = faradaic_current(u_at, V_at)
@printf("starting log-j: V=%.4f, j=%.3e A/m² = %.3f mA/cm²\n", V_at, j0, j0*0.1)

hist_logj = newton_continuation_logj(u_at, V_at, -2.5, build_res, faradaic_current;
                                     n_steps_per_decade = 100,
                                     max_iter = 60, verbose = false)
V_final, _ = hist_logj[end]
@printf("\nlog-j continuation final V = %.4f (%d points)\n", V_final, length(hist_logj))
j_final = faradaic_current(hist_logj[end][2], V_final)
@printf("  j_final = %.3e A/m² = %.3f mA/cm²\n", j_final, j_final*0.1)
