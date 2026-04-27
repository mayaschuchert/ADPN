module Hydrodynamics

# Maps Bloomquist reactor inputs (gap, Q_aq, Q_org) to the Nernst-layer
# thickness δ used by the 1D model and to the Weber numbers used as
# diagnostic flow-regime coordinates. v6 uses a laminar Lévêque correlation
# with no bubble-induced enhancement (deferred to v7).

export ml_min_to_m3_s, d_hydraulic, v_super,
       reynolds, schmidt, sherwood_leveque,
       delta_leveque, weber_numbers,
       W_CHANNEL, L_CHANNEL, NU_KIN, RHO_AQ, RHO_ORG, SIGMA_AN

# ---------- Bloomquist reactor geometry ----------
const W_CHANNEL = 4.0e-3       # m,  FEP gasket channel width (Bloomquist SI)
const L_CHANNEL = 0.16         # m,  serpentine path length = A_active / W_CHANNEL
                               #     (A_active = 6.4 cm², W = 4 mm)

# ---------- Fluid properties (25 °C, water-like) ----------
const NU_KIN    = 1.0e-6       # m² s⁻¹, kinematic viscosity (water at 25 °C)
const RHO_AQ    = 1000.0       # kg m⁻³,  aqueous-phase density
const RHO_ORG   = 810.0        # kg m⁻³,  AN-phase density (Bloomquist SI)
const SIGMA_AN  = 10.5e-3      # N m⁻¹,   AN-water interfacial tension
                               #          (Girifalco-Good, Bloomquist SI)

# Reference D for Schmidt number — AN aqueous diffusivity (the slowest of the
# fitting-relevant neutrals in v6 — see Params.D_aq[6]).
const D_AN_AQ_REF = 2.30e-9    # m² s⁻¹

# ---------- Unit helper ----------
"Convert volumetric flow from mL/min to m³/s."
ml_min_to_m3_s(q_ml_min) = q_ml_min * 1.0e-6 / 60.0

# ---------- Geometry ----------
"Hydraulic diameter d_h = 2·gap·W / (gap + W) for a rectangular channel."
@inline d_hydraulic(gap_m) = 2 * gap_m * W_CHANNEL / (gap_m + W_CHANNEL)

"Superficial velocity v = Q / (gap · W)."
@inline v_super(Q_m3s, gap_m) = Q_m3s / (gap_m * W_CHANNEL)

# ---------- Dimensionless numbers ----------
"Reynolds number Re = v·d_h / ν."
@inline reynolds(v, d_h) = v * d_h / NU_KIN

"Schmidt number Sc = ν / D."
@inline schmidt(D = D_AN_AQ_REF) = NU_KIN / D

"""
    sherwood_leveque(Re, Sc, d_h, L)

Lévêque-developing-BL Sherwood number for laminar flow in a rectangular duct:
Sh = 1.85 · (Re · Sc · d_h / L)^(1/3). Valid for Re·Sc·d_h/L ≫ 1, which holds
across the entire Bloomquist parameter space (Re < 10, Sc ≈ 435).
"""
@inline sherwood_leveque(Re, Sc, d_h, L) = 1.85 * (Re * Sc * d_h / L)^(1/3)

"""
    delta_leveque(gap_m, Q_total_m3s; D_ref = D_AN_AQ_REF) -> δ [m]

Boundary-layer thickness from the Lévêque correlation, using *total* superficial
velocity (sum of phase flows) — see v6 §18.2. Uses AN aqueous diffusivity as
the reference; `D_ref` can be overridden for sensitivity studies but the model
itself does not couple D_AN to δ_lam at runtime.
"""
function delta_leveque(gap_m, Q_total_m3s; D_ref::Float64 = D_AN_AQ_REF)
    d_h = d_hydraulic(gap_m)
    v   = v_super(Q_total_m3s, gap_m)
    Re  = reynolds(v, d_h)
    Sc  = schmidt(D_ref)
    Sh  = sherwood_leveque(Re, Sc, d_h, L_CHANNEL)
    return d_h / Sh
end

"""
    weber_numbers(gap_m, Q_aq_m3s, Q_org_m3s) -> (We_aq, We_org)

Phase Weber numbers We_i = ρ_i · v_i² · gap / σ_AN-water. Diagnostic in v6;
not consumed by the model. Use to overlay Bloomquist points on the regime
map and to cross-check the SI tabulated values to within rounding.
"""
function weber_numbers(gap_m, Q_aq_m3s, Q_org_m3s)
    v_aq  = v_super(Q_aq_m3s,  gap_m)
    v_org = v_super(Q_org_m3s, gap_m)
    We_aq  = RHO_AQ  * v_aq^2  * gap_m / SIGMA_AN
    We_org = RHO_ORG * v_org^2 * gap_m / SIGMA_AN
    return (We_aq = We_aq, We_org = We_org)
end

# v7 stub (left here intentionally as a flagged TODO):
# function delta_actual(gap_m, Q_total_m3s, j_A_m2; bubble_model = nothing)
#     δ_lam = delta_leveque(gap_m, Q_total_m3s)
#     bubble_model === nothing && return δ_lam
#     return δ_lam * f_bubble(j_A_m2, gap_m, Q_total_m3s; bubble_model)
# end

end # module
