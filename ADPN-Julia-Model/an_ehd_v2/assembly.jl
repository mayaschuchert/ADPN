module Assembly

using ..Params
using ..Diffusivity
using ..Chemistry
using ..Kinetics
using ..Transport

export conc_dof, phi_dof, full_residual!, decode_state, faradaic_currents_from_state

# ------------------------------------------------------------
# DOF layout — cell-major, 9 DOFs per cell
#   species k ∈ 1..8  →  idx = 9*(ix-1) + k
#   φ_l               →  idx = 9*ix
# ------------------------------------------------------------
@inline conc_dof(ix::Int, k::Int) = 9 * (ix - 1) + k
@inline phi_dof(ix::Int)          = 9 * ix

# ------------------------------------------------------------
# Decode solution vector -> (c [8×N], phi [N]).
# Parametric on element type so ForwardDiff.Dual numbers propagate through.
# ------------------------------------------------------------
function decode_state(u::AbstractVector{T}, N::Int) where {T<:Real}
    c = zeros(T, 8, N)
    phi = zeros(T, N)
    for ix in 1:N
        base = 9 * (ix - 1)
        @inbounds for k in 1:8
            c[k, ix] = exp(clamp(u[base + k], -50.0, 50.0))
        end
        phi[ix] = u[base + 9]
    end
    return c, phi
end

"""
    faradaic_currents_from_state(u, V, alpha_kin) -> (j1, j2, j3)  [A m⁻²]

Evaluate Tafel currents using the first-cell AN and φ_l as the surface values.
"""
function faradaic_currents_from_state(u::AbstractVector, V::Float64,
                                      alpha_kin::Float64 = 1.0)
    base = 0
    c_AN_surf = exp(clamp(u[base + 6], -50.0, 50.0))
    phi_surf  = u[base + 9]
    return tafel_currents(c_AN_surf, phi_surf, V, alpha_kin)
end

# ------------------------------------------------------------
# Full residual (integrated FV form, §10.5, §10.7)
# ------------------------------------------------------------
"""
    full_residual!(F, u, mesh, eps_org, V, alpha_buf, alpha_kin, c_eq)

Writes the residual F (length 9N) in-place.  Uses integrated form
F[species] = J_L − J_R + S·dx (mol m⁻² s⁻¹) for interior cells, and
current-conservation Σ z_i (J_L − J_R) for φ_l rows.  Bulk cell uses
Dirichlet overrides; electrode face uses Faradaic flux BCs.
"""
function full_residual!(res::AbstractVector{T}, u::AbstractVector{T},
                        mesh, eps_org::Float64, V::Real,
                        alpha_buf::Real, alpha_kin::Real,
                        c_eq) where {T<:Real}

    N = length(mesh.dx)
    @assert length(u) == 9 * N
    @assert length(res) == 9 * N

    # ---- Decode DOFs ----
    c, phi = decode_state(u, N)

    # ---- Interior face fluxes (faces 2..N : between cell ix and ix+1) ----
    J = zeros(T, 8, N + 1)
    @inbounds for ix in 1:(N - 1)
        dx_face = 0.5 * (mesh.dx[ix] + mesh.dx[ix + 1])
        for k in 1:8
            D_k = D_mix(k, eps_org)
            J[k, ix + 1] = sg_flux(c[k, ix], c[k, ix + 1],
                                   phi[ix], phi[ix + 1],
                                   D_k, z_species[k], dx_face)
        end
    end

    # ---- Face 1 (x = 0): Faradaic boundary fluxes ----
    # +x oriented; consumption ⇒ N_i(0) < 0, production ⇒ N_i(0) > 0.
    j1, j2, j3 = tafel_currents(c[6, 1], phi[1], V, alpha_kin)
    J[1, 1] = 0.0                              # H⁺
    J[2, 1] = +(j1 + j2 + j3) / F              # OH⁻ produced
    J[3, 1] = 0.0                              # H₂PO₄⁻
    J[4, 1] = 0.0                              # HPO₄²⁻
    J[5, 1] = 0.0                              # PO₄³⁻
    J[6, 1] = -(2.0 * j1 + j2) / (2.0 * F)     # AN consumed
    J[7, 1] = +(j1) / (2.0 * F)                # ADPN produced
    J[8, 1] = +(j2) / (2.0 * F)                # PN produced
    # Face N+1 is unused (Dirichlet overrides last-cell residual).

    # ---- Assemble residuals ----
    R_buf = zeros(T, 8)
    @inbounds for ix in 1:N
        if ix == N
            # Bulk Dirichlet: force log(c) = log(c_bulk); φ_l = 0 (gauge).
            for k in 1:8
                res[conc_dof(ix, k)] = u[conc_dof(ix, k)] -
                                       log(bulk_concentration(k, c_eq, eps_org))
            end
            res[phi_dof(ix)] = phi[ix]
        else
            buffer_sources!(R_buf,
                            c[1, ix], c[2, ix], c[3, ix], c[4, ix], c[5, ix],
                            alpha_buf)
            dx_ix = mesh.dx[ix]
            for k in 1:8
                res[conc_dof(ix, k)] = J[k, ix] - J[k, ix + 1] + R_buf[k] * dx_ix
            end
            # Current conservation: Σ z_k (J_left − J_right) over charged species (k=1..5).
            s_left  = 0.0
            s_right = 0.0
            for k in 1:5
                s_left  += z_species[k] * J[k, ix]
                s_right += z_species[k] * J[k, ix + 1]
            end
            res[phi_dof(ix)] = s_left - s_right
        end
    end

    return res
end

end # module
