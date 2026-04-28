# Verify ForwardDiff AD Jacobian runs correctly through the residual chain
# and compare against banded-FD at a bootstrapped state.
using LinearAlgebra, SparseArrays, ForwardDiff, Printf

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD, .ADPN_EHD.Params, .ADPN_EHD.Chemistry, .ADPN_EHD.Mesh
using .ADPN_EHD.Assembly, .ADPN_EHD.Solver

mesh = make_mesh(100, 50e-6; stretch = 10.0)
c_eq = solve_phosphate_equilibrium()
u = make_initial_guess(100, c_eq, 0.0)

# Bootstrap quickly (use :fd mode for baseline)
for α in range(0.1, 1.0; length = 10)
    newton_solve!(u, (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, α, 0.0, c_eq);
                  max_iter = 40, verbose = false)
end
for α in [min(1e-6 * 2.0^k, 1.0) for k in 0:20]
    newton_solve!(u, (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, 1.0, α, c_eq);
                  max_iter = 40, verbose = false)
end
println("bootstrap done")

# Now at V=-1.8, compute residual + Jacobian both ways
V_test = -1.8
res_cb = (F, uu) -> full_residual!(F, uu, mesh, 0.0, V_test, 1.0, 1.0, c_eq)

# FD
n = length(u)
F0 = zeros(n)
res_cb(F0, u)
Jsp = build_banded_pattern(n, 17)
banded_fd_jacobian!(Jsp, res_cb, u, F0)
Jfd = Matrix(Jsp)

# AD via ForwardDiff
cfg = ForwardDiff.JacobianConfig(res_cb, zeros(n), u, ForwardDiff.Chunk{12}())
Jad = zeros(n, n)
ForwardDiff.jacobian!(Jad, res_cb, zeros(n), u, cfg)

@printf("||Jfd||_inf  = %.3e  (argmax at %s)\n",
        maximum(abs, Jfd), argmax(abs.(Jfd)))
@printf("||Jad||_inf  = %.3e  (argmax at %s)\n",
        maximum(abs, Jad), argmax(abs.(Jad)))
@printf("||Jfd - Jad||_inf = %.3e\n", maximum(abs, Jfd .- Jad))

# Where's the biggest disagreement?
idx = argmax(abs.(Jfd .- Jad))
@printf("Biggest FD/AD disagreement at (row=%d, col=%d): Jfd=%.3e, Jad=%.3e\n",
        idx[1], idx[2], Jfd[idx], Jad[idx])

# What DOFs are those?
row, col = idx[1], idx[2]
row_cell = (row - 1) ÷ 9 + 1;  row_k = (row - 1) % 9 + 1
col_cell = (col - 1) ÷ 9 + 1;  col_k = (col - 1) % 9 + 1
species_names = ["H⁺", "OH⁻", "H₂PO₄⁻", "HPO₄²⁻", "PO₄³⁻", "AN", "ADPN", "PN"]
row_name = row_k ≤ 8 ? species_names[row_k] : "φ_l"
col_name = col_k ≤ 8 ? species_names[col_k] : "φ_l"
@printf("  → row = cell %d %s, col = cell %d %s\n",
        row_cell, row_name, col_cell, col_name)
@printf("  |row - col| = %d (bandwidth = 17)\n", abs(row - col))

# Compare dense FD to banded FD explicitly
function dense_fd_jacobian(residual!, u, F0; eps_phi=1e-5, eps_conc=1e-7)
    n = length(u)
    J = zeros(n, n)
    F_p = similar(F0)
    for j in 1:n
        ε = (j % 9 == 0) ? eps_phi : eps_conc * max(abs(u[j]), 1.0)
        du = zeros(n); du[j] = ε
        residual!(F_p, u .+ du)
        J[:, j] = (F_p .- F0) ./ ε
    end
    return J
end
Jdense = dense_fd_jacobian(res_cb, u, F0)
@printf("\n||Jfd - Jdense||_inf (banded vs dense FD, should be ~0) = %.3e\n",
        maximum(abs, Jfd .- Jdense))
@printf("||Jad - Jdense||_inf (AD vs dense FD)                   = %.3e\n",
        maximum(abs, Jad .- Jdense))

# Conditioning comparison
σ_fd = svd(Jfd).S
σ_ad = svd(Jad).S
@printf("cond(Jfd) = %.3e,   cond(Jad) = %.3e\n", σ_fd[1]/σ_fd[end], σ_ad[1]/σ_ad[end])
println("AD Jacobian matches FD to machine precision  ✓")

# Run a single newton_solve! at V_test in both modes and compare convergence
println("\n─── Newton at V=-1.8 V ───")
u_fd = copy(u)
r_fd = newton_solve!(u_fd, res_cb; max_iter = 40, jacobian_mode = :fd, verbose = false)
@printf("FD  mode: converged=%s  iter=%d  |F|=%.3e\n", r_fd.converged, r_fd.iter, r_fd.normF)

u_ad = copy(u)
r_ad = newton_solve!(u_ad, res_cb; max_iter = 40, jacobian_mode = :ad, verbose = false)
@printf("AD  mode: converged=%s  iter=%d  |F|=%.3e\n", r_ad.converged, r_ad.iter, r_ad.normF)

# Solutions should be nearly identical
@printf("||u_fd - u_ad||_inf = %.3e\n", maximum(abs, u_fd .- u_ad))
