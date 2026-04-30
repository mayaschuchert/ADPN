module Assembly

using ..Params
using ..Diffusivity
using ..Chemistry
using ..Kinetics
using ..Transport

export conc_dof, phi_dof, full_residual!, decode_state, faradaic_currents_from_state

# Precompute at module load so the type-parameter T in generic functions doesn't shadow
# the Params.T temperature constant.
const _RT_F = R_gas * T / F     # thermal voltage [V] ≈ 0.02569 V

# ------------------------------------------------------------
# DOF layout — cell-major, 10 DOFs per cell
#   species k ∈ 1..9  →  idx = 10*(ix-1) + k
#   φ_l               →  idx = 10*ix
# ------------------------------------------------------------
@inline conc_dof(ix::Int, k::Int) = 10 * (ix - 1) + k
@inline phi_dof(ix::Int)          = 10 * ix

# ------------------------------------------------------------
# Decode solution vector -> (c [9×N], phi [N]).
# Parametric on element type so ForwardDiff.Dual numbers propagate through.
# ------------------------------------------------------------
function decode_state(u::AbstractVector{T}, N::Int) where {T<:Real}
    c = zeros(T, 9, N)
    phi = zeros(T, N)
    for ix in 1:N
        base = 10 * (ix - 1)
        @inbounds for k in 1:9
            c[k, ix] = exp(clamp(u[base + k], -50.0, 50.0))
        end
        phi[ix] = u[base + 10]
    end
    return c, phi
end

"""
    faradaic_currents_from_state(u, V, alpha_kin) -> (j1, j2, j3, j4)  [A m⁻²]

Evaluate Tafel currents using the first-cell AN and φ_l as the surface values.
"""
function faradaic_currents_from_state(u::AbstractVector, V::Float64,
                                      alpha_kin::Float64 = 1.0)
    base = 0
    c_AN_surf = exp(clamp(u[base + 6], -50.0, 50.0))
    phi_surf  = u[base + 10]
    return tafel_currents(c_AN_surf, phi_surf, V, alpha_kin)
end

# ------------------------------------------------------------
# Full residual (integrated FV form, §10.5, §10.7)
# ------------------------------------------------------------
"""
    full_residual!(F, u, mesh, eps_org, V, alpha_buf, alpha_kin, c_eq)

Writes the residual F (length 10N) in-place.  Uses integrated form
F[species] = J_L − J_R + S·dx (mol m⁻² s⁻¹) for interior cells, and
current-conservation Σ z_i (J_L − J_R) for φ_l rows.  Bulk cell uses
Dirichlet overrides; electrode face uses Faradaic flux BCs.

For Lees-style S1.11 sequential fitting, partial currents at the electrode
face can be clamped to experimental values via the existing `Kinetics.FIXED_J`
Ref (set NaN → BV-computed; finite → use directly). The substitution point is
inside `tafel_currents`, so this residual function needs no separate kwarg.
"""
function full_residual!(res::AbstractVector{T}, u::AbstractVector{T},
                        mesh, eps_org::Float64, V::Real,
                        alpha_buf::Real, alpha_kin::Real,
                        c_eq) where {T<:Real}

    N = length(mesh.dx)
    @assert length(u) == 10 * N
    @assert length(res) == 10 * N

    # ---- Decode DOFs ----
    c, phi = decode_state(u, N)

    # ---- Interior face fluxes (faces 2..N : between cell ix and ix+1) ----
    J    = zeros(T, 9, N + 1)
    J_Na = zeros(T, N + 1)    # Na⁺ flux; J_Na[1]=0 (no electrode reaction for Na⁺)
    D_Na_mix = D_Na_aq * (1.0 - eps_org)
    @inbounds for ix in 1:(N - 1)
        dx_face = 0.5 * (mesh.dx[ix] + mesh.dx[ix + 1])
        for k in 1:9
            D_k = D_mix(k, eps_org)
            J[k, ix + 1] = sg_flux(c[k, ix], c[k, ix + 1],
                                   phi[ix], phi[ix + 1],
                                   D_k, z_species[k], dx_face)
        end
        # Na⁺ concentration inferred from electroneutrality at each cell
        c_Na_L = max(c[2, ix]   + c[3, ix]   + 2.0 * c[4, ix]   + 3.0 * c[5, ix]   - c[1, ix],   1e-30)
        c_Na_R = max(c[2,ix+1]  + c[3,ix+1]  + 2.0 * c[4,ix+1]  + 3.0 * c[5,ix+1]  - c[1,ix+1],  1e-30)
        J_Na[ix + 1] = sg_flux(c_Na_L, c_Na_R, phi[ix], phi[ix + 1], D_Na_mix, 1, dx_face)
    end

    # ---- Face 1 (x = 0): Faradaic boundary fluxes ----
    # +x oriented; consumption ⇒ N_i(0) < 0, production ⇒ N_i(0) > 0.
    # Per-branch FIXED_J override (NaN → BV; finite → clamp) is handled inside
    # tafel_currents — no extra logic needed here.
    j1, j2, j3, j4 = tafel_currents(c[6, 1], phi[1], V, alpha_kin)
    J[1, 1] = 0.0                                    # H⁺
    J[2, 1] = +(j1 + j2 + j3 + j4) / F              # OH⁻ produced
    J[3, 1] = 0.0                                    # H₂PO₄⁻
    J[4, 1] = 0.0                                    # HPO₄²⁻
    J[5, 1] = 0.0                                    # PO₄³⁻
    J[6, 1] = -(2.0 * j1 + j2 + 3.0 * j4) / (2.0 * F)  # AN consumed (3 AN per TCH)
    J[7, 1] = +(j1) / (2.0 * F)                     # ADPN produced
    J[8, 1] = +(j2) / (2.0 * F)                     # PN produced
    J[9, 1] = +(j4) / (2.0 * F)                     # TCH produced
    # Face N+1 is unused (Dirichlet overrides last-cell residual).

    # ---- Assemble residuals ----
    R_buf = zeros(T, 9)
    @inbounds for ix in 1:N
        if ix == N
            # Bulk Dirichlet: force log(c) = log(c_bulk); φ_l = 0 (gauge).
            for k in 1:9
                res[conc_dof(ix, k)] = u[conc_dof(ix, k)] -
                                       log(bulk_concentration(k, c_eq, eps_org))
            end
            res[phi_dof(ix)] = phi[ix]
        else
            buffer_sources!(R_buf,
                            c[1, ix], c[2, ix], c[3, ix], c[4, ix], c[5, ix],
                            alpha_buf)
            dx_ix = mesh.dx[ix]
            for k in 1:9
                res[conc_dof(ix, k)] = J[k, ix] - J[k, ix + 1] + R_buf[k] * dx_ix
            end
            # φ equation: total current conservation including Na⁺ background cation.
            # In steady state the divergence of total ionic current is zero (no charge source).
            # Σ z_k J_k is constant across all faces; residual = left_face − right_face = 0.
            s_left  = J_Na[ix]
            s_right = J_Na[ix + 1]
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
