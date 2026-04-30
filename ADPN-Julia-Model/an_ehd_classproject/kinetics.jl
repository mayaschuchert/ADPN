module Kinetics

using ..Params

export tafel_currents,
       set_kinetic_override!, clear_kinetic_override!, with_kinetic_override,
       FIXED_J

# Optional override of (j₀,r, α_c,r, n_r) used by the fit workflow. When
# `nothing` (the default), `tafel_currents` reads from Params and uses default
# reaction orders; Stage 1/2/2m/3 paths are unaffected. Stage 4v3 sets the
# Ref before each fixed-j solve so the optimiser can vary kinetics without
# touching Params.
#
# v7: extended to 4 reactions (added TCH). j0/ac are 4-tuples; n is a 3-tuple
# (n_ADN, n_PN, n_TCH). HER (j0[3], ac[3]) is frozen in Stage 4v3 — the
# override still carries the full 4-tuple so solve_at_j has a uniform interface.
# Per-row sequential-fit override. NaN sentinel = compute from Tafel; finite = use directly.
# Fixed values bypass alpha_kin scaling (always applied at full strength).
const FIXED_J = Ref{NTuple{4,Float64}}((NaN, NaN, NaN, NaN))

const KIN_OVERRIDE = Ref{Union{Nothing,
                               @NamedTuple{j0::NTuple{4,Float64},
                                           ac::NTuple{4,Float64},
                                           n::NTuple{3,Float64}}}}(nothing)

set_kinetic_override!(j0::NTuple{4,Float64}, ac::NTuple{4,Float64},
                      n::NTuple{3,Float64} = (2.0, 1.0, 3.0)) =
    (KIN_OVERRIDE[] = (j0 = j0, ac = ac, n = n); nothing)
clear_kinetic_override!() = (KIN_OVERRIDE[] = nothing; nothing)

function with_kinetic_override(f, j0::NTuple{4,Float64}, ac::NTuple{4,Float64},
                               n::NTuple{3,Float64} = (2.0, 1.0, 3.0))
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
        -> (j1, j2, j3, j4)

Cathodic Tafel current densities [A m⁻²] for ADPN, PN, HER, TCH at the electrode.
Sign convention: positive = cathodic. Overpotential η_r = (φ_s − φ_l) − E⁰_r.

j₀_r are in SI (A m⁻²). If the v7 fit override is active, j₀_r and α_c,r come
from the override tuple instead of Params constants. Default: read Params.
"""
@inline function tafel_currents(c_AN_surface::Real,
                                phi_l_surface::Real,
                                phi_s::Real,
                                alpha_kin::Real = 1.0)

    ov = KIN_OVERRIDE[]
    if ov === nothing
        j01, j02, j03, j04 = j0_1, j0_2, j0_3, j0_4
        a1, a2, a3, a4      = alpha_c1, alpha_c2, alpha_c3, alpha_c4
        n_ADN, n_PN, n_TCH  = 2.0, 1.0, n_TCH_default
    else
        j01, j02, j03, j04 = ov.j0
        a1, a2, a3, a4      = ov.ac
        n_ADN, n_PN, n_TCH  = ov.n
    end

    eta1 = (phi_s - phi_l_surface) - E0_1
    eta2 = (phi_s - phi_l_surface) - E0_2
    eta3 = (phi_s - phi_l_surface) - E0_3
    eta4 = (phi_s - phi_l_surface) - E0_TCH

    # Guard against negative surface AN (numerical noise)
    cA = max(c_AN_surface, zero(c_AN_surface))

    j1 = j01 * (cA / c_ref)^n_ADN * exp(-a1 * F * eta1 / (R_gas * T))
    j2 = j02 * (cA / c_ref)^n_PN  * exp(-a2 * F * eta2 / (R_gas * T))
    j3 = j03 *                      exp(-a3 * F * eta3 / (R_gas * T))
    j4 = j04 * (cA / c_ref)^n_TCH * exp(-a4 * F * eta4 / (R_gas * T))

    fj = FIXED_J[]
    return (isnan(fj[1]) ? alpha_kin * j1 : fj[1],
            isnan(fj[2]) ? alpha_kin * j2 : fj[2],
            isnan(fj[3]) ? alpha_kin * j3 : fj[3],
            isnan(fj[4]) ? alpha_kin * j4 : fj[4])
end

end # module
