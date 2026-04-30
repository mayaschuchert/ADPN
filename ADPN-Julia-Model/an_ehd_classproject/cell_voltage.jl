module CellVoltage

# v6 §17 — external cell-voltage decomposition.
#
#   V_cell  =  V_CE  +  |V_cathode_SHE|  +  j · R_series        (positive magnitude)
#   R_series = (gap − δ) / κ_eff(c_bulk, ε_org)  +  R_contact
#   κ_eff    = κ_dilute(c_bulk) · (1 − ε_org)^1.5               (Bruggeman)
#   κ_dilute = (F²/RT) · Σ_i z_i² · D_i,aq · c_i,bulk            (Newman §11.3)
#
# v6 freezes (V_CE, R_contact) at literature defaults — see guide §20.5.
# Bubble void Bruggeman correction (1 − ε_gas)^1.5 is deferred to v7.

using ..Params

export V_CE_DEFAULT, R_CONTACT_DEFAULT,
       D_NA_AQ,
       kappa_dilute, kappa_eff,
       R_series, V_cell_predicted, V_cathode_target

# ---------- v6 defaults for the frozen voltage scalars ----------
const V_CE_DEFAULT       = 1.7      # V vs SHE  — lumped E°_OER + ⟨η_anode⟩ on SS
const R_CONTACT_DEFAULT  = 1.0e-4   # Ω·m²      — spring-probe + Cd-foil contact stack

# Na⁺ aqueous diffusivity (Params.D_aq is indexed 1..8 over the 8 transported
# species; Na⁺ is recovered from electroneutrality and not stored there). v6
# §17.2 uses the CRC value for κ_dilute.
const D_NA_AQ = 1.33e-9             # m² s⁻¹  (CRC Handbook, 25 °C)

# ---------- Dilute-solution conductivity ----------
"""
    kappa_dilute(c_eq) -> κ [S m⁻¹]

Dilute-solution conductivity from the bulk equilibrium tuple returned by
`Chemistry.solve_phosphate_equilibrium`. Uses Newman-style sum
κ = (F²/RT)·Σ z_i²·D_i,aq·c_i with ions {H⁺, OH⁻, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻, Na⁺}.

`c_eq` must have fields `H, OH, H2PO4, HPO4, PO4, Na` in mol m⁻³.
"""
function kappa_dilute(c_eq)
    coeff = F^2 / (R_gas * T)
    # Indices in Params.D_aq: 1=H⁺, 2=OH⁻, 3=H₂PO₄⁻, 4=HPO₄²⁻, 5=PO₄³⁻
    return coeff * (
        1.0 * D_aq[1] * c_eq.H     +    # z=+1
        1.0 * D_aq[2] * c_eq.OH    +    # z=-1
        1.0 * D_aq[3] * c_eq.H2PO4 +    # z=-1
        4.0 * D_aq[4] * c_eq.HPO4  +    # z=-2 → z² = 4
        9.0 * D_aq[5] * c_eq.PO4   +    # z=-3 → z² = 9
        1.0 * D_NA_AQ * c_eq.Na          # z=+1, Na⁺ from electroneutrality
    )
end

"""
    kappa_eff(c_eq, eps_org) -> κ_eff [S m⁻¹]

Effective conductivity of the two-phase (or single-phase) electrolyte after
the Bruggeman porosity correction for organic-droplet volume fraction:
κ_eff = κ_dilute · (1 − ε_org)^1.5. Reduces to κ_dilute below ε_sat where
v6's regime-aware D_mix already turns off; the Bruggeman factor accounts for
ionic-conduction blockage by organic droplets in the two-phase regime.
"""
@inline kappa_eff(c_eq, eps_org) = kappa_dilute(c_eq) * (1.0 - eps_org)^1.5

# ---------- Series resistance ----------
"""
    R_series(gap_m, delta_m, c_eq, eps_org; R_contact=R_CONTACT_DEFAULT) -> Ω·m²

Bulk-electrolyte ohmic resistance over (gap − δ) plus contact resistance.
The diffusion-layer ohmic drop (0 ≤ x ≤ δ on the cathode side) is *already*
in V_cathode_SHE via the model's φ_l(0) — do not double-count.
"""
function R_series(gap_m::Real, delta_m::Real, c_eq, eps_org::Real;
                  R_contact::Float64 = R_CONTACT_DEFAULT)
    @assert gap_m > delta_m "δ must be smaller than the gap; got δ = $delta_m, gap = $gap_m"
    return (gap_m - delta_m) / kappa_eff(c_eq, eps_org) + R_contact
end

# ---------- V_cell forward / inverse ----------
"""
    V_cell_predicted(V_cathode_SHE, j_A_m2, gap_m, delta_m, eps_org, c_eq;
                     V_CE=V_CE_DEFAULT, R_contact=R_CONTACT_DEFAULT) -> V [V]

Positive cell-voltage magnitude (matches Bloomquist's convention).
`V_cathode_SHE` is the model's solved cathode potential (negative).
`j_A_m2` is total cathodic current density [A m⁻²] (use `j_mA_cm2 * 10`).
"""
function V_cell_predicted(V_cathode_SHE::Real, j_A_m2::Real,
                          gap_m::Real, delta_m::Real,
                          eps_org::Real, c_eq;
                          V_CE::Float64       = V_CE_DEFAULT,
                          R_contact::Float64  = R_CONTACT_DEFAULT)
    Rs = R_series(gap_m, delta_m, c_eq, eps_org; R_contact = R_contact)
    return V_CE + abs(V_cathode_SHE) + j_A_m2 * Rs
end

"""
    V_cathode_target(V_cell_meas, j_A_m2, gap_m, delta_m, eps_org, c_eq;
                     V_CE=V_CE_DEFAULT, R_contact=R_CONTACT_DEFAULT) -> V [V vs SHE]

Inverse of `V_cell_predicted`: given a measured `V_cell` magnitude, return the
(signed, negative) cathode potential the model should solve at. Useful for
forward-prediction sanity checks; v6's actual fit uses fixed-j (not fixed-V),
so this is diagnostic.
"""
function V_cathode_target(V_cell_meas::Real, j_A_m2::Real,
                          gap_m::Real, delta_m::Real,
                          eps_org::Real, c_eq;
                          V_CE::Float64       = V_CE_DEFAULT,
                          R_contact::Float64  = R_CONTACT_DEFAULT)
    Rs = R_series(gap_m, delta_m, c_eq, eps_org; R_contact = R_contact)
    return -(V_cell_meas - V_CE - j_A_m2 * Rs)
end

end # module
