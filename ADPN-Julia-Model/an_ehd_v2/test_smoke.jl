# Smoke tests — verify module wiring & basic identities before running stage 1.
using Printf

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.Params
using .ADPN_EHD.Chemistry
using .ADPN_EHD.Mesh
using .ADPN_EHD.Transport
using .ADPN_EHD.Diffusivity
using .ADPN_EHD.Assembly

println("── Smoke tests ──")

# 1. solve_phosphate_equilibrium
c_eq = solve_phosphate_equilibrium()
@printf("equilibrium: pH = %.4f, cH = %.3e, cOH = %.3e, cHPO4 = %.3e, cPO4 = %.3e\n",
        c_eq.pH, c_eq.H, c_eq.OH, c_eq.HPO4, c_eq.PO4)

@assert abs(c_eq.H * c_eq.OH - K_w) / K_w < 1e-8                 "K_w mismatch"
@assert abs(c_eq.H * c_eq.HPO4 / c_eq.H2PO4 - K_a2) / K_a2 < 1e-8 "K_a2 mismatch"
@assert abs(c_eq.H * c_eq.PO4  / c_eq.HPO4  - K_a3) / K_a3 < 1e-8 "K_a3 mismatch"

# 2. Buffer residual at equilibrium
R = zeros(9)
buffer_sources!(R, c_eq.H, c_eq.OH, c_eq.H2PO4, c_eq.HPO4, c_eq.PO4)
@printf("max|R_buf| at equilibrium = %.3e\n", maximum(abs.(R)))
@assert maximum(abs.(R)) < 1e-6 "buffer residual at equilibrium too large"

# 3. Mesh
mesh = make_mesh(100, 50e-6; stretch = 10.0)
@printf("mesh: dx_min = %.3e, dx_max = %.3e, ratio = %.3f, Σdx = %.3e (should = 5e-5)\n",
        minimum(mesh.dx), maximum(mesh.dx), maximum(mesh.dx)/minimum(mesh.dx), sum(mesh.dx))
@assert abs(sum(mesh.dx) - 50e-6) / 50e-6 < 1e-12

# 4. SG flux neutral species
F_neutral = sg_flux(10.0, 5.0, 0.0, 0.0, 1e-9, 0, 1e-6)
F_ref     = -1e-9 * (5.0 - 10.0) / 1e-6
@printf("SG neutral: %.6e  (expected %.6e)\n", F_neutral, F_ref)
@assert abs(F_neutral - F_ref) < 1e-20

# 5. SG flux at small potential → approaches centered diff
F_ion_zero = sg_flux(10.0, 5.0, 0.0, 0.0, 1e-9, -1, 1e-6)
@printf("SG ion (Δφ=0): %.6e  (expected %.6e)\n", F_ion_zero, F_ref)
@assert abs(F_ion_zero - F_ref) / abs(F_ref) < 1e-8

# 6. D_mix (regime-aware)
@printf("EPS_ORG_SAT = %.4f  (expected ≈ 0.0862, Convention A)\n", EPS_ORG_SAT)
@printf("D_mix(AN, 0.02, single-phase)  = %.3e   (expected D_aq = 2.30e-9)\n", D_mix(6, 0.02))
@printf("D_mix(AN, 0.05, single-phase)  = %.3e   (expected D_aq = 2.30e-9)\n", D_mix(6, 0.05))
@printf("D_mix(AN, 0.15, two-phase)     = %.3e   (expected 2.86e-9)\n", D_mix(6, 0.15))
@printf("D_mix(AN, 0.30, two-phase)     = %.3e   (expected 3.41e-9)\n", D_mix(6, 0.30))
@printf("D_mix(OH, 0.30, two-phase)     = %.3e   (expected 3.69e-9)\n", D_mix(2, 0.30))
@assert abs(D_mix(6, 0.02) - 2.30e-9) < 1e-12
@assert abs(D_mix(6, 0.30) - (0.3 * 6.0e-9 + 0.7 * 2.3e-9)) < 1e-12
@assert abs(D_mix(2, 0.30) - 0.7 * 5.27e-9) < 1e-12

# 6b. c_AN_bulk regime sweep
@printf("\n-- c_AN_bulk vs ε_org --\n")
for eps in (0.02, 0.05, 0.08, 0.09, 0.15, 0.25, 0.30)
    regime = eps < EPS_ORG_SAT ? "single-phase" : "two-phase"
    @printf("ε_org = %.2f   c_AN,bulk = %7.1f mol/m³   (%s)\n",
            eps, c_AN_bulk(eps), regime)
end
@assert c_AN_bulk(0.02) < 400.0
@assert c_AN_bulk(0.15) ≈ C_AN_SAT
@assert c_AN_bulk(0.30) ≈ C_AN_SAT

# 7. Initial guess ⇒ zero residual with α_buf = α_kin = 0 at a PHYSICAL ε_org.
# (ε_org = 0 now gives c_AN,bulk = 0 which is truly no-AN; a well-defined
#  reference state needs ε_org > 0. Use the Stage 1 anchor 0.02.)
eps_test = 0.02
u0 = make_initial_guess(100, c_eq, eps_test)
F0 = zeros(length(u0))
full_residual!(F0, u0, mesh, eps_test, -1.0, 0.0, 0.0, c_eq)
@printf("|F|∞ at α_buf=α_kin=0, ε_org=%.2f: %.3e (should be ~0)\n",
        eps_test, maximum(abs.(F0)))
@assert maximum(abs.(F0)) < 1e-10 "zero-state residual too large"

# 8. At α_buf=1, α_kin=0, still should be near zero (equilibrium IC is exact buffer eq)
fill!(F0, 0.0)
full_residual!(F0, u0, mesh, eps_test, -1.0, 1.0, 0.0, c_eq)
@printf("|F|∞ at α_buf=1, α_kin=0, ε_org=%.2f: %.3e (should be ~0)\n",
        eps_test, maximum(abs.(F0)))

println("── smoke tests passed ──")
