module Kinetics

using ..Params

export tafel_currents,
       set_kinetic_override!, clear_kinetic_override!, with_kinetic_override

# Optional override of (j₀,r, α_c,r, n_r) used by the v6 fit workflow. When
# `nothing` (the default), `tafel_currents` reads j₀,r and α_c,r from Params
# and uses the historical fixed AN reaction orders n_1=2, n_2=1; all Stage
# 1 / 2 / 2m / 3 paths are unaffected. Stage 4 sets the Ref before each
# fixed-j solve so the optimiser can vary kinetics without touching Params.
#
# v6.x: `n` was added to the override to promote the AN reaction orders
# (n_1 for ADPN, n_2 for PN) to fit parameters. Backwards-compat: callers
# that still pass only (j0, ac) get the historical defaults n=(2.0, 1.0).
const KIN_OVERRIDE = Ref{Union{Nothing,
                               @NamedTuple{j0::NTuple{3,Float64},
                                           ac::NTuple{3,Float64},
                                           n::NTuple{2,Float64}}}}(nothing)

set_kinetic_override!(j0::NTuple{3,Float64}, ac::NTuple{3,Float64},
                      n::NTuple{2,Float64} = (2.0, 1.0)) =
    (KIN_OVERRIDE[] = (j0 = j0, ac = ac, n = n); nothing)
clear_kinetic_override!() = (KIN_OVERRIDE[] = nothing; nothing)

function with_kinetic_override(f, j0::NTuple{3,Float64}, ac::NTuple{3,Float64},
                               n::NTuple{2,Float64} = (2.0, 1.0))
    prev = KIN_OVERRIDE[]
    set_kinetic_override!(j0, ac, n)
    try
        return f()
    finally
        KIN_OVERRIDE[] = prev
    end
end

"""
    tafel_currents(c_AN_surface, phi_l_surface, phi_s, alpha_kin=1.0)
        -> (j1, j2, j3)

Cathodic Tafel current densities [A m⁻²] for ADPN, PN, HER at the electrode.
Sign convention: positive = cathodic. Overpotential η_r = (φ_s − φ_l) − E⁰_r.

j₀_r are in SI (A m⁻²), already converted from mA cm⁻² in Params. If the
v6 fit override is active (`KIN_OVERRIDE[] !== nothing`), j₀_r and α_c,r come
from the override tuple instead of Params constants. Default: read Params.
"""
@inline function tafel_currents(c_AN_surface::Real,
                                phi_l_surface::Real,
                                phi_s::Real,
                                alpha_kin::Real = 1.0)

    ov = KIN_OVERRIDE[]
    if ov === nothing
        j01, j02, j03 = j0_1, j0_2, j0_3
        a1, a2, a3    = alpha_c1, alpha_c2, alpha_c3
        n1, n2        = 2.0, 1.0
    else
        j01, j02, j03 = ov.j0
        a1, a2, a3    = ov.ac
        n1, n2        = ov.n
    end

    eta1 = (phi_s - phi_l_surface) - E0_1
    eta2 = (phi_s - phi_l_surface) - E0_2
    eta3 = (phi_s - phi_l_surface) - E0_3

    # Guard against negative surface AN (numerical noise)
    cA = max(c_AN_surface, zero(c_AN_surface))

    j1 = j01 * (cA / c_ref)^n1 * exp(-a1 * F * eta1 / (R_gas * T))
    j2 = j02 * (cA / c_ref)^n2 * exp(-a2 * F * eta2 / (R_gas * T))
    j3 = j03 *                   exp(-a3 * F * eta3 / (R_gas * T))

    return (alpha_kin * j1, alpha_kin * j2, alpha_kin * j3)
end

end # module
