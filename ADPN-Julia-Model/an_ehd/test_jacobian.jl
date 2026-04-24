# Diagnostic: compare banded-FD Jacobian to a dense-FD Jacobian at a state where
# continuation fails. We reproduce the bootstrap state at V = -1.0, then try the
# first V step and inspect the Jacobian.
using LinearAlgebra, SparseArrays, Printf

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.Params
using .ADPN_EHD.Chemistry
using .ADPN_EHD.Mesh
using .ADPN_EHD.Assembly
using .ADPN_EHD.Solver

# dense FD Jacobian (slow, diagnostic only)
function dense_fd_jacobian(residual!, u, F0; eps_phi=1e-5, eps_conc=1e-7)
    n = length(u)
    J = zeros(n, n)
    F_pert = similar(F0)
    for j in 1:n
        ε = (j % 9 == 0) ? eps_phi : eps_conc * max(abs(u[j]), 1.0)
        du = zeros(n); du[j] = ε
        residual!(F_pert, u .+ du)
        J[:, j] = (F_pert .- F0) ./ ε
    end
    return J
end

mesh = make_mesh(100, 50e-6; stretch = 10.0)
c_eq = solve_phosphate_equilibrium()
u = make_initial_guess(100, c_eq, 0.0)

# Replay the bootstrap
println("─ Bootstrap ─")
for α in range(0.1, 1.0; length = 10)
    res_cb = (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, α, 0.0, c_eq)
    r = newton_solve!(u, res_cb)
    # (result ignored)
end
α_list = [min(1e-6 * 2.0^k, 1.0) for k in 0:20]; α_list[end] = 1.0
for α in α_list
    res_cb = (F, u) -> full_residual!(F, u, mesh, 0.0, -1.0, 1.0, α, c_eq)
    r = newton_solve!(u, res_cb)
end
println("bootstrap done")

# Now go to V = -1.675 (last successful step in the live run) and try -1.7
function diag_loop(u)
for V in [-1.0, -1.2, -1.4, -1.6, -1.65, -1.675, -1.7, -1.75, -1.8]
    res_cb = (F, uu) -> full_residual!(F, uu, mesh, 0.0, V, 1.0, 1.0, c_eq)
    F0 = zeros(length(u))
    res_cb(F0, u)
    normF = maximum(abs.(F0))

    # Try to continue — but snapshot state BEFORE Newton mutates it.
    u_snap = copy(u)
    u_try  = copy(u)
    r = newton_solve!(u_try, res_cb; verbose=false, max_iter=25)
    if r.converged
        u = u_try
        @printf("V=%+.3f: starting |F|=%.3e → converged in %d iter, |F|=%.3e\n",
                V, normF, r.iter, r.normF)
    else
        @printf("V=%+.3f: starting |F|=%.3e → FAILED, iter=%d, |F|=%.3e\n",
                V, normF, r.iter, r.normF)

        # Diagnostic: compare banded vs dense Jacobian AT INITIAL state.
        F0_snap = zeros(length(u_snap))
        res_cb(F0_snap, u_snap)
        Jsp = build_banded_pattern(length(u_snap), 17)
        banded_fd_jacobian!(Jsp, res_cb, u_snap, F0_snap)
        Jdense = dense_fd_jacobian(res_cb, u_snap, F0_snap)

        # Difference
        Jbnd_dense = Matrix(Jsp)
        diff_max = maximum(abs.(Jbnd_dense .- Jdense))
        @printf("  max |Jbanded - Jdense| = %.3e\n", diff_max)
        # Find rows/cols where they differ most
        for _ in 1:1
            idx = argmax(abs.(Jbnd_dense .- Jdense))
            @printf("  max diff at (row %d, col %d): Jbnd=%.3e Jdense=%.3e\n",
                    idx[1], idx[2], Jbnd_dense[idx], Jdense[idx])
        end

        # rank of dense
        σ = svd(Jdense).S
        @printf("  Jdense singular values: σ_min=%.3e, σ_max=%.3e, cond=%.3e\n",
                minimum(σ), maximum(σ), maximum(σ)/minimum(σ))
        # report rank
        rnk_tol = 1e-10 * maximum(σ)
        r_rank  = count(>(rnk_tol), σ)
        @printf("  dense rank (tol=%.3e): %d / %d\n", rnk_tol, r_rank, length(σ))

        # try dense solve
        try
            du = -(Jdense \ F0)
            @printf("  dense solve OK: max|du|=%.3e\n", maximum(abs.(du)))
        catch e
            @printf("  dense solve also fails: %s\n", sprint(showerror, e))
        end
        break
    end
end
end

diag_loop(u)
