module Params

export F, R_gas, T, K_w, K_a2, K_a3,
       k1f, k2f, k3f, k1r, k2r, k3r,
       D_aq, D_org, m_partition, z_species, n_species,
       C_AN_SAT, c_ref,
       M_AN, RHO_AN, m_AN, MOLAR_DENSITY_AN,
       E0_1, E0_2, E0_3,
       alpha_c1, alpha_c2, alpha_c3,
       j0_1, j0_2, j0_3,
       C_P_total, c_Na_input,
       EPS_ORG_SAT

# ---------- Physical constants ----------
const F      = 96_485.332      # C mol⁻¹
const R_gas  = 8.314463        # J mol⁻¹ K⁻¹
const T      = 298.15          # K (25 °C)

# ---------- Equilibrium constants (SI: mol m⁻³) ----------
# All K values here are in mol/m³. Literature M-based values:
#   K_w  = 10⁻¹⁴ M² × (1000 mol/m³ per M)² = 10⁻⁸ (mol/m³)²
#   Ka2  =  6.3×10⁻⁸ M × 1000                =  6.3×10⁻⁵ mol/m³
#   Ka3  =  4.5×10⁻¹³ M × 1000               =  4.5×10⁻¹⁰ mol/m³
# Ka2 and Ka3 are pinned by the OH-pathway K_eq values below
# (K_eq,R2 = Ka2/K_w = 6.3×10⁶ M⁻¹; K_eq,R3 = Ka3/K_w = 45 M⁻¹).
const K_w  = 1.0e-8            # (mol/m³)², pKw = 14.00 in M-units
const K_a2 = 6.3e-5            # mol/m³, pKa2 = 7.20
const K_a3 = 4.5e-10           # mol/m³, pKa3 = 12.35

# ---------- Buffer rate constants (OH⁻-pathway formulation) ----------
# R1:  H₂O  ⇌  H⁺ + OH⁻
#       k₁,f = 1.4×10⁻³ M/s    = 1.4 mol/(m³·s)          (Eigen–De Maeyer 1955)
#       k₁,r = k₁,f / K_w      (M⁻¹s⁻¹ → m³/(mol·s))
#
# R2:  OH⁻ + H₂PO₄⁻  ⇌  HPO₄²⁻ + H₂O
#       k₂,f = 1×10⁸  M⁻¹ s⁻¹  = 1×10⁵ m³/(mol·s)
#       K_eq,R2 = Ka2/K_w = 6.3×10⁶ M⁻¹ = 6.3×10³ m³/mol
#       k₂,r    = k₂,f / K_eq,R2                          (units: s⁻¹)
#
# R3:  OH⁻ + HPO₄²⁻  ⇌  PO₄³⁻  + H₂O
#       k₃,f = 2×10⁶  M⁻¹ s⁻¹  = 2×10³ m³/(mol·s)
#       K_eq,R3 = Ka3/K_w = 45 M⁻¹ = 4.5×10⁻² m³/mol
#       k₃,r    = k₃,f / K_eq,R3                          (units: s⁻¹)

const k1f = 1.4                # mol/(m³·s)             — zeroth-order in solutes
const k1r = k1f / K_w          # m³/(mol·s)    = 1.4×10⁸

# OH⁻-pathway constants (R2 and R3 are second-order forward, first-order reverse)
const K_eq_R2 = K_a2 / K_w     # m³/mol        = 6.3×10³
const K_eq_R3 = K_a3 / K_w     # m³/mol        = 4.5×10⁻²

const k2f = 1.0e5              # m³/(mol·s)
const k3f = 2.0e3              # m³/(mol·s)
const k2r = k2f / K_eq_R2      # s⁻¹           ≈ 15.87
const k3r = k3f / K_eq_R3      # s⁻¹           ≈ 4.44×10⁴

# ---------- Species ordering ----------
# 1:H⁺  2:OH⁻  3:H₂PO₄⁻  4:HPO₄²⁻  5:PO₄³⁻  6:AN  7:ADPN  8:PN
const n_species = 8
const z_species = (+1, -1, -1, -2, -3, 0, 0, 0)

# ---------- Diffusivities (m² s⁻¹, 25 °C) ----------
const D_aq = [9.31e-9, 5.27e-9, 0.846e-9, 0.690e-9, 0.610e-9,
              2.30e-9, 1.50e-9, 2.30e-9]
const D_org = [0.0, 0.0, 0.0, 0.0, 0.0,
               6.00e-9, 3.90e-9, 6.00e-9]

# Partition coefficients m_i = c_i,org / c_i,aq at local equilibrium.
# Ions are insoluble in the organic phase (m_i = 0). Neutrals from
# Suwanvaipattana et al. 2017; ADPN & PN not directly tabulated — use
# dinitrile analogue (SN, 7.72) for ADPN and AN-analog (11.59) for PN.
# Used ONLY by the m_i-corrected D_eff formulation (§4.2).
#                    H+   OH-  H2PO4- HPO4²- PO4³-   AN     ADPN   PN
const m_partition = [0.0, 0.0, 0.0,   0.0,   0.0,    11.59, 7.72,  11.59]

# ---------- AN solubility / partition ----------
# ε_org is the volume fraction of AN *added* per total volume — a loading
# parameter, not a droplet volume fraction. In the single-phase regime all
# AN dissolves in the aqueous phase; in the two-phase regime the aqueous
# phase saturates and excess AN forms organic droplets.
const M_AN              = 0.05306       # kg/mol        (acrylonitrile, C₃H₃N)
const RHO_AN            = 806.0         # kg/m³          (neat AN, 25 °C)
const m_AN              = 11.59         # —              (partition coefficient, Suwanvaipattana 2017)
const MOLAR_DENSITY_AN  = RHO_AN / M_AN # ≈ 15 191 mol/m³ (molar density of neat AN)
const C_AN_SAT          = MOLAR_DENSITY_AN / m_AN  # ≈ 1 310 mol/m³ (aqueous saturation)
# Single-phase ⇄ two-phase threshold, Convention A (per total solution volume):
# ε_org · ρ_AN/M_AN = C_AN_SAT
const EPS_ORG_SAT       = C_AN_SAT / MOLAR_DENSITY_AN                # ≈ 0.0862

const c_ref       = 1000.0     # mol m⁻³ (= 1 M) for Tafel

# ---------- Standard potentials (V vs SHE) ----------
const E0_1 = -1.30
const E0_2 = -1.30
const E0_3 = -0.83

# ---------- Transfer coefficients ----------
const alpha_c1 = 0.5
const alpha_c2 = 0.5
const alpha_c3 = 0.4

# ---------- Exchange current densities (A m⁻²) ----------
# Guide gives mA cm⁻²; convert ×10.
const j0_1 = 1.0e-4 * 10.0     # 10⁻³ A m⁻²
const j0_2 = 1.0e-4 * 10.0     # 10⁻³ A m⁻²
const j0_3 = 1.0e-6 * 10.0     # 10⁻⁵ A m⁻²

# ---------- Bulk electrolyte composition (Modestino group) ----------
const C_P_total  = 500.0       # mol m⁻³ total phosphate (0.5 M Na₃PO₄)
const c_Na_input = 1520.0      # 3×500 (Na⁺) + 20 (TBA⁺) lumped, mol m⁻³

end # module
