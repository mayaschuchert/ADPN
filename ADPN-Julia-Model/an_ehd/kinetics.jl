module Kinetics

using ..Params

export tafel_currents

"""
    tafel_currents(c_AN_surface, phi_l_surface, phi_s, alpha_kin=1.0)
        -> (j1, j2, j3)

Cathodic Tafel current densities [A m⁻²] for ADPN, PN, HER at the electrode.
Sign convention: positive = cathodic. Overpotential η_r = (φ_s − φ_l) − E⁰_r.

j₀_r are in SI (A m⁻²), already converted from mA cm⁻² in Params.
"""
@inline function tafel_currents(c_AN_surface::Real,
                                phi_l_surface::Real,
                                phi_s::Real,
                                alpha_kin::Real = 1.0)

    eta1 = (phi_s - phi_l_surface) - E0_1
    eta2 = (phi_s - phi_l_surface) - E0_2
    eta3 = (phi_s - phi_l_surface) - E0_3

    # Guard against negative surface AN (numerical noise)
    cA = max(c_AN_surface, zero(c_AN_surface))

    j1 = j0_1 * (cA / c_ref)^2 * exp(-alpha_c1 * F * eta1 / (R_gas * T))
    j2 = j0_2 * (cA / c_ref)    * exp(-alpha_c2 * F * eta2 / (R_gas * T))
    j3 = j0_3 *                   exp(-alpha_c3 * F * eta3 / (R_gas * T))

    return (alpha_kin * j1, alpha_kin * j2, alpha_kin * j3)
end

end # module
