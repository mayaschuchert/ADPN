# Acrylonitrile Electrohydrodimerization (EHD) — 1D Planar Electrode Model

**Bui Lab | NYU Tandon School of Engineering | April 2026 | Guide v5**

Scope: 1D Nernst diffusion layer (δ), Nernst–Planck transport with migration, Tafel kinetics for ADPN/PN/HER, phosphate buffer chemistry (OH⁻-pathway), regime-aware multiphase effective diffusivity on a Cd cathode. **Operating temperature: T = 298.15 K (25 °C).**

---

## Convention

> All concentrations c_i in this guide are **aqueous-phase** values [mol m⁻³ of aqueous solution]. Under the local-equilibrium assumption (Da >> 1), organic and aqueous phases equilibrate instantaneously at every position: c_i,org(x) = m_i × c_i,aq(x). The organic phase acts as a parallel transport pathway captured entirely through the effective diffusivity D_i,mix. **No explicit phase transfer term R_PT appears** in the governing equations. No volume-fraction prefactors on any source terms.

> **ε_org convention.** ε_org is the volume fraction of AN added per unit total solution volume — a *loading* parameter, not strictly a droplet volume fraction. Below the solubility threshold ε_sat it represents AN fully dissolved in a single aqueous phase; above ε_sat it represents AN in excess of the saturation limit, forming organic droplets. **AN concentrations use Convention A — moles per total solution volume** (see §7.2).

---

## Table of Contents

1. Physical Domain
2. Species and Degrees of Freedom
3. Governing Equations
4. Mixture-Averaged Diffusivities (regime-aware)
5. Electrochemical Kinetics
6. Phosphate Buffer Chemistry (OH⁻-pathway)
7. Boundary Conditions
8. Three-Parameter Sweep
9. Parameter Tables
10. Numerical Methods
11. Solution Caching
12. Implementation Stages
13. Physicality Checks
14. Module Structure
15. Pitfalls
16. Potential Referencing
17. Experimental Data
18. Fitting Strategy
19. Required Plots

---

## 1. Physical Domain

The model solves steady-state species transport across a stagnant Nernst diffusion layer of thickness δ [m] adjacent to a planar cadmium cathode (at x = 0), with the well-mixed bulk electrolyte at x = δ. Framework: Corpus et al. (Joule 2023) and Weng, Bell & Weber (PCCP 2018).

The diffusion layer contains dispersed organic droplets characterised by their volume fraction ε_org [—]. Under the local-equilibrium assumption, the organic and aqueous phases are in instantaneous equilibrium everywhere — the droplets provide a parallel diffusion pathway for species that partition into the organic phase. All electrochemical reactions occur at the electrode surface and enter the model as flux boundary conditions at x = 0. H₂ desorbs immediately and is not tracked.

```
 ELECTRODE (Cd)         DIFFUSION LAYER (δ)                   BULK ELECTROLYTE
     x = 0          organic droplets in local equil.               x = δ
                    |                                    |
 Tafel kinetics     |  NP transport with migration       |  Dirichlet BCs
 ADPN + PN + HER    |  D_i,mix(ε_org) eff. diffusivity  |  c_i = c_i,bulk
 Flux BCs           |  Phosphate buffer (OH-pathway)     |
                    |  No R_PT (local equilibrium)       |
                    |  Current conservation (φ_l)        |
```

Key findings from Bloomquist et al. (CEJ 2026, 528, 172125):

- FE_ADPN increases from <50% to >80% when ε_org exceeds the solubility limit (~0.086)
- Bubble-induced convection dominates over inlet flow regime
- High selectivity maintained at j > 200 mA cm⁻²
- ε_org explored from 0.02 to 0.30

---

## 2. Species and Degrees of Freedom

The model tracks 9 dissolved species. Eight are independent unknowns; Na⁺ is computed from electroneutrality.

| Index | Species | Charge z_i | Role |
|-------|---------|------------|------|
| 1 | H⁺ | +1 | Proton |
| 2 | OH⁻ | −1 | Hydroxide (product of all 3 cathodic reactions) |
| 3 | H₂PO₄⁻ | −1 | Buffer |
| 4 | HPO₄²⁻ | −2 | Buffer |
| 5 | PO₄³⁻ | −3 | Buffer |
| 6 | AN (CH₂=CHCN) | 0 | Acrylonitrile (reactant) |
| 7 | ADPN (NC(CH₂)₄CN) | 0 | Adiponitrile (product) |
| 8 | PN (CH₃CH₂CN) | 0 | Propionitrile (byproduct) |
| 9 | Na⁺ | +1 | From electroneutrality (not a DOF) |

> **TBA⁺ treatment:** Tetrabutylammonium (0.02 M = 20 mol m⁻³) is present as a selectivity agent but is lumped into the Na⁺ background — both carry z = +1, and since Na⁺ is recovered from electroneutrality rather than transported, no mobility mismatch enters the model. The effective background Na⁺ input is c_Na,input = 3 × 500 + 20 = **1520 mol m⁻³** (see §6.4).

Na⁺ is determined algebraically from local electroneutrality:

```
c_Na = c_OH + c_H₂PO₄ + 2·c_HPO₄ + 3·c_PO₄ − c_H
```

Each finite-volume cell carries **9 degrees of freedom**: 8 log-concentrations (species 1–8) plus the electrolyte potential φ_l [V]. Log-concentrations are defined as u_i = ln(c_i) and recovered via c_i = exp(clamp(u_i, −50, 50)). With N_mesh = 100: **900 total unknowns**.

---

## 3. Governing Equations

### 3.1 Nernst–Planck Transport with Migration

For each species i at steady state:

$$0 = \frac{\partial}{\partial x}\!\left[D_{i,\mathrm{mix}}\frac{\partial c_i}{\partial x} + \frac{z_i F}{RT}D_{i,\mathrm{mix}}\,c_i\frac{\partial\phi_\ell}{\partial x}\right] + R_{\mathrm{buf},i}$$

```
0 = d/dx [D_i,mix × dc_i/dx + (z_i·F)/(R·T) × D_i,mix × c_i × dφ_l/dx] + R_buf,i
```

| Symbol | Definition | Unit |
|--------|-----------|------|
| c_i | Aqueous-phase concentration of species i | mol m⁻³ |
| D_i,mix | Effective diffusivity through emulsion (§4, regime-aware) | m² s⁻¹ |
| z_i | Charge number | — |
| F | Faraday constant = 96,485.332 | C mol⁻¹ |
| R | Gas constant = 8.314463 | J mol⁻¹ K⁻¹ |
| T | Temperature = 298.15 K (25 °C) | K |
| φ_l | Electrolyte potential (solved DOF) | V |
| R_buf,i | Buffer reaction source for species i (§6) | mol m⁻³ s⁻¹ |
| x | Spatial coordinate (0 = electrode, δ = bulk) | m |
| δ | Diffusion layer thickness | m |

**No R_PT term.** Under the local-equilibrium assumption, the organic and aqueous phases equilibrate instantaneously and there is no driving force for interphase transfer. All multiphase effects are captured through D_i,mix.

**No prefactor on R_buf.** Buffer reactions occur in the aqueous phase, the rate law uses aqueous concentrations, and c_i is an aqueous concentration. No correction is needed.

### 3.2 Current Conservation (Potential Equation)

The 9th equation per cell enforces zero net ionic current in the diffusion layer. No external current flows through the solution between x = 0 and x = δ — current enters and exits only at the electrode. This constraint determines the electrolyte potential φ_l at every position:

$$\frac{\partial}{\partial x}\!\left[\sum_i z_i N_i\right] = 0$$

where N_i is the Nernst–Planck flux of species i and the sum runs over charged species only. In the FV discretisation, this becomes: the sum of ionic SG fluxes at the right face of each cell must equal the sum at the left face.

> **Well-posedness (gauge fix):** The current-conservation equation `d/dx[Σ z_i N_i] = 0` is invariant under a uniform shift of φ_l — the system has one null direction. The Dirichlet condition φ_l(δ) = 0 (§7.2) kills this null direction and makes the Jacobian non-singular. Without it, Newton cannot converge.

### 3.3 Scharfetter–Gummel Flux Discretisation

Standard centred-difference flux approximations become numerically unstable when the dimensionless migration Péclet number α = z_i F Δφ/(RT) exceeds ~2 per cell, producing spurious oscillations. The SG scheme uses an exponential fitting to handle arbitrarily large α stably.

$$J_i = -\frac{D_{i,\mathrm{mix}}}{\Delta x}\left[B(\alpha)\,c_R - B(-\alpha)\,c_L\right], \quad \alpha = \frac{z_i F(\phi_R - \phi_L)}{RT}$$

with the Bernoulli functions:

$$B(\alpha) = \frac{\alpha}{e^\alpha - 1}, \qquad B(-\alpha) = \frac{\alpha\,e^\alpha}{e^\alpha - 1}$$

where c_L, c_R are aqueous concentrations in the left and right cells [mol m⁻³]; φ_L, φ_R are electrolyte potentials [V]; and Δx is the distance between cell centres [m]. α is clamped to [−700, 700] to avoid floating-point overflow.

> **Taylor expansion for small α (required for AD).** At α = 0 exactly (e.g. the equilibrium initial guess), the naïve Bernoulli form divides by zero. A hard `if abs(α) < 1e-10: return centered_diff` branch *also* fails for automatic differentiation (AD) — it creates a spurious zero derivative at the α = 0 seed point because the centered-difference form has no φ dependence. **Use the Taylor series for |α| < 0.01** to keep the function smooth *and* its φ-derivative continuous through α = 0:
>
> ```
> B(α)  = 1 − α/2 + α²/12 − α⁴/720 + ...    (for |α| < 0.01)
> B(−α) = 1 + α/2 + α²/12 + α⁴/720 + ...
> ```

```julia
function sg_flux(c_L, c_R, phi_L, phi_R, D_mix, z_i, dx, T=298.15)
    z_i == 0 && return -D_mix * (c_R - c_L) / dx
    alpha = clamp(z_i * F * (phi_R - phi_L) / (R * T), -700.0, 700.0)
    if abs(alpha) < 0.01
        a2 = alpha * alpha
        B_pos = 1.0 - alpha / 2.0 + a2 / 12.0
        B_neg = 1.0 + alpha / 2.0 + a2 / 12.0
    else
        ea  = exp(alpha)
        den = ea - 1.0
        B_pos = alpha / den
        B_neg = alpha * ea / den
    end
    return -D_mix / dx * (B_pos * c_R - B_neg * c_L)
end
```

---

## 4. Mixture-Averaged Diffusivities — Regime-Aware

### 4.1 Two regimes

AN loading partitions the electrolyte into two regimes:

| Regime | Threshold | Physical state | D_i,mix |
|---|---|---|---|
| **Single-phase** | ε_org < ε_sat | All AN dissolved in aqueous phase — no droplets | **D_i,mix = D_i,aq** (no mixing) |
| **Two-phase** | ε_org ≥ ε_sat | Aqueous phase saturated at C_AN_SAT; excess AN forms organic droplets | **D_i,mix = ε_org · D_i,org + (1 − ε_org) · D_i,aq** (arithmetic mean) |

The threshold ε_sat is determined by the aqueous saturation concentration of AN (§7.2):

$$\varepsilon_{\mathrm{sat}} = \frac{C_{\mathrm{AN,SAT}}}{\rho_{\mathrm{AN}}/M_{\mathrm{AN}}} = \frac{1310}{15{,}191} \approx \mathbf{0.0862}$$

Below ε_sat there are no organic droplets at all; the "organic diffusivity pathway" is physically absent. Only above ε_sat does the volume-averaged diffusivity become meaningful.

Ions are insoluble in the organic phase, so D_i,org = 0 for all ions. The two-phase formula becomes:

$$D_{i,\mathrm{mix}}^{\,\mathrm{ion}} = (1 - \varepsilon_{\mathrm{org}})\,D_{i,\mathrm{aq}} \quad (\varepsilon_{\mathrm{org}} \ge \varepsilon_{\mathrm{sat}})$$

Increasing ε_org above ε_sat **reduces ionic transport** (ohmic penalty) while **enhancing neutral-species transport** (AN parallel pathway). Neutral organic diffusivities are estimated via the Stokes–Einstein viscosity ratio: D_i,org ≈ D_i,aq × (μ_aq/μ_org) ≈ 2.6 × D_i,aq (μ_H₂O = 0.89 mPa·s, μ_AN = 0.34 mPa·s at 25°C).

| Species | D_aq (m²/s) | D_org (m²/s) | Notes |
|---------|-------------|--------------|-------|
| AN | 2.30 × 10⁻⁹ | 6.00 × 10⁻⁹ | Stokes–Einstein ×2.6 |
| ADPN | 1.50 × 10⁻⁹ | 3.90 × 10⁻⁹ | Stokes–Einstein ×2.6 |
| PN | 2.30 × 10⁻⁹ | 6.00 × 10⁻⁹ | Stokes–Einstein ×2.6 |
| H⁺ | 9.31 × 10⁻⁹ | 0 | Insoluble in organic |
| OH⁻ | 5.27 × 10⁻⁹ | 0 | Insoluble in organic |
| H₂PO₄⁻ | 0.846 × 10⁻⁹ | 0 | Insoluble in organic |
| HPO₄²⁻ | 0.690 × 10⁻⁹ | 0 | Insoluble in organic |
| PO₄³⁻ | 0.610 × 10⁻⁹ | 0 | Insoluble in organic |
| Na⁺ | 1.33 × 10⁻⁹ | 0 | Insoluble in organic |

D_aq from Suwanvaipattana et al. (2017) and CRC Handbook. All values at T = 25 °C.

```julia
const T   = 298.15   # K, operating temperature

const D_aq = [9.31e-9, 5.27e-9, 0.846e-9, 0.69e-9, 0.61e-9,   # H⁺ OH⁻ H₂PO₄⁻ HPO₄²⁻ PO₄³⁻
              2.3e-9, 1.5e-9, 2.3e-9]                            # AN ADPN PN
const D_org = [0.0, 0.0, 0.0, 0.0, 0.0,                         # ions: insoluble
               6.0e-9, 3.9e-9, 6.0e-9]                           # AN ADPN PN

const EPS_ORG_SAT = C_AN_SAT / MOLAR_DENSITY_AN    # ≈ 0.0862

function D_mix(i, eps_org)
    if eps_org < EPS_ORG_SAT
        return D_aq[i]                                  # single-phase: no organic
    else
        return eps_org * D_org[i] + (1.0 - eps_org) * D_aq[i]   # two-phase: arithmetic
    end
end
```

### Effect of ε_org on Mixture Diffusivities

| ε_org | Regime | D_AN,mix (m²/s) | D_OH,mix (m²/s) |
|-------|--------|----------------|----------------|
| 0.02 | single-phase | 2.30 × 10⁻⁹ | 5.27 × 10⁻⁹ |
| 0.05 | single-phase | 2.30 × 10⁻⁹ | 5.27 × 10⁻⁹ |
| 0.08 | single-phase | 2.30 × 10⁻⁹ | 5.27 × 10⁻⁹ |
| **0.0862** | **threshold** | **2.62 × 10⁻⁹** | **4.82 × 10⁻⁹** |
| 0.09 | two-phase | 2.63 × 10⁻⁹ | 4.80 × 10⁻⁹ |
| 0.15 | two-phase | 2.86 × 10⁻⁹ | 4.48 × 10⁻⁹ |
| 0.25 | two-phase | 3.23 × 10⁻⁹ | 3.95 × 10⁻⁹ |
| 0.30 | two-phase | 3.41 × 10⁻⁹ | 3.69 × 10⁻⁹ |

The transition at ε_sat is a ~14% step in D_AN (jumping from pure D_aq to the arithmetic mean). This reflects the physical onset of organic droplets providing a parallel transport pathway.

AN transport improves by 48% at ε_org = 0.30 (relative to single-phase D_aq) while OH⁻ transport degrades by 30% — this quantifies the selectivity–conductivity trade-off in the two-phase regime.

### 4.2 m_i-Corrected Upgrade (Future)

If the arithmetic mean is too weak to reproduce the experimental FE enhancement, the physically correct local-equilibrium effective diffusivity includes the partition coefficient m_i = c_i,org/c_i,aq:

$$D_{i,\mathrm{eff}} = (1-\varepsilon_{\mathrm{org}})\,D_{i,\mathrm{aq}} + \varepsilon_{\mathrm{org}}\,m_i\,D_{i,\mathrm{org}} \quad (\varepsilon_{\mathrm{org}} \ge \varepsilon_{\mathrm{sat}})$$

For AN (m_AN = 11.59) at ε_org = 0.25: D_AN,eff = 19.1 × 10⁻⁹ m²/s vs D_AN,mix = 3.23 × 10⁻⁹. For ions: m_i = 0, so D_ion,eff = (1−ε_org) D_i,aq — no change.

```julia
# One-line switch from arithmetic to m_i-corrected (two-phase only):
D_mix(i, eo) = eo < EPS_ORG_SAT ? D_aq[i] : eo * D_org[i] + (1-eo) * D_aq[i]           # arithmetic
D_eff(i, eo) = eo < EPS_ORG_SAT ? D_aq[i] : eo * m[i] * D_org[i] + (1-eo) * D_aq[i]   # m-corrected
```

---

## 5. Electrochemical Kinetics

### 5.1 Competing Cathodic Reactions

| Reaction | Stoichiometry | n_e |
|----------|---------------|-----|
| R1: Electrohydrodimerization | 2 AN + 2 H₂O + 2e⁻ → ADPN + 2 OH⁻ | 2 |
| R2: Hydrogenation | AN + 2 H₂O + 2e⁻ → PN + 2 OH⁻ | 2 |
| R3: Hydrogen evolution | 2 H₂O + 2e⁻ → H₂ + 2 OH⁻ | 2 |

All three reactions transfer n_e = 2 electrons and produce 2 OH⁻.

### 5.2 Tafel Rate Expressions

All reactions operate far from equilibrium; the cathodic Tafel approximation applies. The current densities j_r [A m⁻²] at the electrode surface are:

**R1 — ADPN** (second-order in AN):

$$j_1 = j_{0,1}\!\left(\frac{c_{\mathrm{AN}}}{c_{\mathrm{ref}}}\right)^{2}\exp\!\left(-\frac{\alpha_{c,1} F\,\eta_1}{RT}\right)$$

**R2 — PN** (first-order in AN):

$$j_2 = j_{0,2}\,\frac{c_{\mathrm{AN}}}{c_{\mathrm{ref}}}\exp\!\left(-\frac{\alpha_{c,2} F\,\eta_2}{RT}\right)$$

**R3 — HER** (independent of AN):

$$j_3 = j_{0,3}\exp\!\left(-\frac{\alpha_{c,3} F\,\eta_3}{RT}\right)$$

The overpotential for each reaction: η_r = (φ_s - φ_l) - E⁰_r.

| Symbol | Definition | Unit |
|--------|-----------|------|
| j_r | Current density for reaction r (positive = cathodic) | A m⁻² |
| j₀,r | Exchange current density for reaction r | A m⁻² |
| c_AN | Aqueous AN concentration at the electrode surface | mol m⁻³ |
| c_ref | Reference concentration = 1,000 mol m⁻³ (= 1 M) | mol m⁻³ |
| α_c,r | Cathodic transfer coefficient for reaction r | — |
| η_r | Overpotential for reaction r | V |
| φ_s | Applied cathode potential | V vs SHE |
| φ_l | Solved electrolyte potential at the electrode | V |
| E⁰_r | Standard reduction potential for reaction r | V vs SHE |

> **Unit note:** j₀,r is given in the fitting table (§9.2) in mA cm⁻². Convert to SI via ×10: 1 mA cm⁻² = 10 A m⁻². c_ref = 1000 mol m⁻³ so that (c_AN/c_ref) is dimensionless; c_AN can reach C_AN_SAT ≈ 1310 mol m⁻³ at saturation, which is physical.

The c² vs c dependence is the central selectivity feature: as c_AN drops at the electrode under high current density, ADPN selectivity falls faster than PN.

The Faradaic efficiency for ADPN:

$$\mathrm{FE}_{\mathrm{ADPN}} = \frac{j_1}{j_1 + j_2 + j_3} \times 100\%$$

---

## 6. Phosphate Buffer Chemistry — OH⁻-Pathway

### 6.1 Equilibria

Three equilibria are tracked. H₃PO₄ dissociation (pK = 2.12) is omitted because it is negligible at the operating pH of 12–13. Water autoprotolysis uses the classical Eigen–De Maeyer rate; phosphate deprotonations use the **OH⁻-pathway** forms, which are the physically dominant kinetic paths at high pH.

| # | Reaction (forward direction) | K_eq (25°C) |
|---|---|---|
| R1 | H₂O ⇌ H⁺ + OH⁻ | K_w = 10⁻¹⁴ M² = 10⁻⁸ (mol/m³)² |
| R2 | **OH⁻ + H₂PO₄⁻ ⇌ HPO₄²⁻ + H₂O** | K_eq,R2 = Ka2/K_w = 6.3 × 10⁶ M⁻¹ = 6.3 × 10³ m³/mol |
| R3 | **OH⁻ + HPO₄²⁻ ⇌ PO₄³⁻ + H₂O** | K_eq,R3 = Ka3/K_w = 45 M⁻¹ = 4.5 × 10⁻² m³/mol |

Equivalent acid-dissociation constants (consistent with the OH⁻-pathway K_eq values):

| Symbol | Value | Unit |
|---|---|---|
| K_w | 10⁻⁸ | (mol/m³)² |
| Ka2 | K_eq,R2 · K_w = 6.3 × 10⁻⁵ | mol/m³ (pKa2 = 7.20) |
| Ka3 | K_eq,R3 · K_w = 4.5 × 10⁻¹⁰ | mol/m³ (pKa3 = 12.35) |

> **Why OH⁻-pathway, not H⁺-pathway?** At pH ≈ 13, c_H ≈ 10⁻¹⁰ mol/m³ and c_OH ≈ 10² mol/m³. The classical form `H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻` requires a reverse rate `k_r · c_H · c_HPO₄` that offsets the forward rate through a tiny c_H — driving k_r to unphysically huge values (~10¹⁰ m³/mol/s). The OH⁻-pathway `OH⁻ + H₂PO₄⁻ → HPO₄²⁻ + H₂O` keeps c_OH (large) in the forward rate and gives physically reasonable rate constants. Both give identical thermodynamic equilibrium.

Reverse rate constants are derived from K_eq: **k_r = k_f / K_eq** (never set independently).

### 6.2 Rate Expressions

For each equilibrium, the net rate [mol m⁻³ s⁻¹] is forward minus reverse:

$$r_1 = k_{1,f} \;-\; k_{1,r}\,c_{\mathrm{H}^+}\,c_{\mathrm{OH}^-}$$

$$r_2 = k_{2,f}\,c_{\mathrm{OH}^-}\,c_{\mathrm{H_2PO_4^-}} \;-\; k_{2,r}\,c_{\mathrm{HPO_4^{2-}}}$$

$$r_3 = k_{3,f}\,c_{\mathrm{OH}^-}\,c_{\mathrm{HPO_4^{2-}}} \;-\; k_{3,r}\,c_{\mathrm{PO_4^{3-}}}$$

**Species allocation (OH⁻-pathway stoichiometry):**

$$R_{\mathrm{buf}}[\mathrm{H}^+] = r_1 \qquad R_{\mathrm{buf}}[\mathrm{OH}^-] = r_1 - r_2 - r_3$$
$$R_{\mathrm{buf}}[\mathrm{H_2PO_4^-}] = -r_2 \qquad R_{\mathrm{buf}}[\mathrm{HPO_4^{2-}}] = r_2 - r_3 \qquad R_{\mathrm{buf}}[\mathrm{PO_4^{3-}}] = r_3$$

R_buf = 0 for AN, ADPN, and PN. **No prefactor** — c_i is aqueous; R_buf is rate per aqueous volume. Charge conservation Σ z_i · S_i = 0 holds for each individual reaction.

```julia
function buffer_sources!(R, c_H, c_OH, c_H2PO4, c_HPO4, c_PO4, alpha_buf=1.0)
    r1 = k1f - k1r * c_H * c_OH                             # H₂O ⇌ H⁺ + OH⁻
    r2 = k2f * c_OH * c_H2PO4 - k2r * c_HPO4                # OH⁻ + H₂PO₄⁻ ⇌ HPO₄²⁻ + H₂O
    r3 = k3f * c_OH * c_HPO4  - k3r * c_PO4                 # OH⁻ + HPO₄²⁻ ⇌ PO₄³⁻  + H₂O
    R[1] = alpha_buf * r1                     # H⁺     — water only
    R[2] = alpha_buf * (r1 - r2 - r3)         # OH⁻    — water + consumed by phosphates
    R[3] = alpha_buf * (-r2)                   # H₂PO₄⁻
    R[4] = alpha_buf * (r2 - r3)              # HPO₄²⁻
    R[5] = alpha_buf * r3                     # PO₄³⁻
    R[6:8] .= 0.0                              # AN, ADPN, PN
end
```

The optional `alpha_buf` argument (default 1.0) is the buffer continuation parameter used during bootstrap (see §10.6, §12).

### 6.3 Rate Constants (SI units)

| Constant | Value | Unit | Source / Derivation |
|---|---|---|---|
| k₁,f | **1.4** | mol/(m³·s) | Eigen–De Maeyer 1955; 1.4 × 10⁻³ M/s × 1000 |
| k₁,r | **1.4 × 10⁸** | m³/(mol·s) | = k₁,f / K_w (≡ 1.4 × 10¹¹ M⁻¹ s⁻¹) |
| k₂,f | **1 × 10⁵** | m³/(mol·s) | 10⁸ M⁻¹ s⁻¹ / 1000 (second-order forward) |
| k₂,r | **≈ 15.87** | s⁻¹ | = k₂,f / K_eq,R2 (unimolecular reverse) |
| k₃,f | **2 × 10³** | m³/(mol·s) | 2 × 10⁶ M⁻¹ s⁻¹ / 1000 |
| k₃,r | **≈ 4.44 × 10⁴** | s⁻¹ | = k₃,f / K_eq,R3 |

Dimensional consistency check: `mol/(m³·s) − m³/(mol·s) · (mol/m³) · (mol/m³)` = `mol/(m³·s)` for R1; `m³/(mol·s) · (mol/m³)² − s⁻¹ · (mol/m³)` = `mol/(m³·s)` for R2, R3. ✓

### 6.4 Bulk Equilibrium Initial Guess

**This is the most critical step for convergence.** Buffer reactions are stiff: forward and reverse rates are both large but nearly cancel at equilibrium. An initial guess that is even slightly off-equilibrium produces buffer source terms of O(10³–10⁵) mol m⁻³ s⁻¹ — orders of magnitude larger than transport fluxes — causing Newton to fail or stall even at α_buf = 0.001.

**The fix: solve the phosphate charge balance exactly before starting Newton.**

The equilibrium speciation follows from two constraints:

**Phosphate conservation:**
$$c_{\mathrm{H_2PO_4}} + c_{\mathrm{HPO_4}} + c_{\mathrm{PO_4}} = C_{P,\mathrm{total}}$$

Expressing each species in terms of c_H alone:
$$c_{\mathrm{HPO_4}} = \frac{C_{P,\mathrm{total}}}{c_H/K_{a2} + 1 + K_{a3}/c_H}, \quad c_{\mathrm{H_2PO_4}} = \frac{c_H \, c_{\mathrm{HPO_4}}}{K_{a2}}, \quad c_{\mathrm{PO_4}} = \frac{K_{a3} \, c_{\mathrm{HPO_4}}}{c_H}$$

**Charge balance** (one equation in c_H):
$$c_{\mathrm{Na,input}} + c_H = \frac{K_w}{c_H} + c_{\mathrm{H_2PO_4}} + 2\,c_{\mathrm{HPO_4}} + 3\,c_{\mathrm{PO_4}}$$

where c_Na,input = 3 × C_Na₃PO₄ + C_TBA-OH = 3 × 500 + 20 = **1520 mol m⁻³** (TBA⁺ lumped into Na⁺, see §2).

Solve via 1D bisection on log₁₀(c_H) over pH ∈ [6, 15]:

```julia
const K_w  = 1.0e-8    # (mol/m³)², water autoprotolysis at 25°C
const K_a2 = 6.3e-5    # mol/m³, H₂PO₄⁻ ↔ H⁺ + HPO₄²⁻ at 25°C
const K_a3 = 4.5e-10   # mol/m³, HPO₄²⁻ ↔ H⁺ + PO₄³⁻ at 25°C

function solve_phosphate_equilibrium(C_P_total=500.0, c_Na_input=1520.0)
    function charge_balance(log10_cH)
        cH    = 10.0^log10_cH
        cOH   = K_w / cH
        denom = cH / K_a2 + 1.0 + K_a3 / cH
        cHPO4  = C_P_total / denom
        cH2PO4 = cH * cHPO4 / K_a2
        cPO4   = K_a3 * cHPO4 / cH
        return c_Na_input + cH - (cOH + cH2PO4 + 2.0*cHPO4 + 3.0*cPO4)
    end
    log10_cH_eq = bisect(charge_balance, -12.0, -3.0)   # or Roots.find_zero
    cH = 10.0^log10_cH_eq
    ...
    return (H=cH, OH=cOH, H2PO4=cH2PO4, HPO4=cHPO4, PO4=cPO4,
            Na=c_Na_input, pH=-log10(cH/1000.0))
end
```

**Verification — always run before starting Newton:**

```julia
c_eq = solve_phosphate_equilibrium()
@assert abs(c_eq.H * c_eq.OH              - K_w)  / K_w  < 1e-8   "K_w mismatch"
@assert abs(c_eq.H * c_eq.HPO4 / c_eq.H2PO4 - K_a2) / K_a2 < 1e-8   "Ka2 mismatch"
@assert abs(c_eq.H * c_eq.PO4  / c_eq.HPO4  - K_a3) / K_a3 < 1e-8   "Ka3 mismatch"

R_check = zeros(8)
buffer_sources!(R_check, c_eq.H, c_eq.OH, c_eq.H2PO4, c_eq.HPO4, c_eq.PO4)
@assert maximum(abs.(R_check)) < 1e-6  "Buffer not at equilibrium"
```

**Expected equilibrium values** for 0.5 M Na₃PO₄ + 0.02 M TBA-OH at 25°C with OH⁻-pathway constants:

| Species | Value | Unit |
|---------|---|------|
| c_H | ≈ 9.4 × 10⁻¹¹ | mol m⁻³ |
| c_OH | ≈ 106 | mol m⁻³ |
| c_H₂PO₄ | ≈ 1.3 × 10⁻⁴ | mol m⁻³ |
| c_HPO₄ | ≈ 86 | mol m⁻³ |
| c_PO₄ | ≈ 414 | mol m⁻³ |
| c_Na | 1520 | mol m⁻³ |
| pH | ≈ 13.03 | — |

Most phosphate sits as PO₄³⁻ since pH ≈ 13.03 is well above pKa3 = 12.35.

**Building the full DOF initial guess vector:**

```julia
function make_initial_guess(N_mesh, c_eq, eps_org; c_seed=1e-3)
    u0 = zeros(9 * N_mesh)
    log_c0 = [log(c_eq.H),    log(c_eq.OH),
              log(c_eq.H2PO4), log(c_eq.HPO4), log(c_eq.PO4),
              log(c_AN_bulk(eps_org)),      # see §7.2 — requires ε_org > 0
              log(c_seed),                   # ADPN: small seed
              log(c_seed)]                   # PN:   small seed
    for ix in 1:N_mesh
        u0[9*(ix-1)+1 : 9*(ix-1)+8] .= log_c0
        u0[9*(ix-1)+9] = 0.0
    end
    return u0
end
```

With this initialisation and α_buf = 0 (buffer off), α_kin = 0 (kinetics off), the residual is **identically zero**: no concentration gradients, no current, no buffer sources. Newton converges in 0 iterations — the ideal starting point for the bootstrap protocol (§12). Requires ε_org > 0 so that c_AN,bulk > 0 (the ε_org = 0 limit is now truly pathological — see §7.2).

---

## 7. Boundary Conditions

### 7.1 Electrode (x = 0): Faradaic Flux

The species fluxes at the electrode surface are set by Faraday's law applied to the three cathodic reactions. Two AN molecules are consumed per ADPN produced and one AN per PN produced:

$$N_{\mathrm{AN}}\big|_{x=0} = -\frac{2j_1 + j_2}{2F}$$
$$N_{\mathrm{ADPN}}\big|_{x=0} = +\frac{j_1}{2F} \qquad N_{\mathrm{PN}}\big|_{x=0} = +\frac{j_2}{2F}$$
$$N_{\mathrm{OH}^-}\big|_{x=0} = +\frac{j_1 + j_2 + j_3}{F}$$

All other species are non-electroactive: N_i(0) = 0 for H⁺ and phosphates.

The electrolyte potential φ_l(0) is **not prescribed** — it is determined by the current conservation equation (§3.2). The overpotentials in the Tafel rate expressions use this solved value.

### 7.2 Bulk (x = δ): Dirichlet — AN via **Convention A**

At the bulk boundary, all species are fixed at their well-mixed bulk values and the electrolyte potential is referenced to zero:

$$c_i(\delta) = c_{i,\mathrm{bulk}} \qquad \phi_\ell(\delta) = 0$$

The bulk equilibrium concentrations of H⁺, OH⁻, and the three phosphate species are computed by `solve_phosphate_equilibrium()` (§6.4). ADPN and PN bulk concentrations are set to small seed values (~10⁻³ mol m⁻³) to avoid log(0) in the initial guess.

**AN bulk concentration (Convention A — moles per total solution volume):**

AN is loaded at a volume fraction ε_org relative to the total solution volume. With ρ_AN = 806 kg/m³ (neat AN density), M_AN = 53.06 g/mol, and m_AN = 11.59 (partition coefficient), the molar density of neat AN is ρ_AN/M_AN = 15,191 mol/m³, and the aqueous saturation concentration is C_AN_SAT = ρ_AN/(M_AN·m_AN) = **1,310 mol/m³**.

| Regime | Condition | c_AN,bulk | Physical picture |
|---|---|---|---|
| **Single-phase** | ε_org < ε_sat ≈ 0.0862 | **ε_org · ρ_AN/M_AN** | All AN dissolved; c_AN linearly rises with loading |
| **Two-phase** | ε_org ≥ ε_sat | **C_AN_SAT = 1310 mol/m³** | Aqueous saturated; excess AN in droplets |

Saturation threshold: ε_sat · (ρ_AN/M_AN) = C_AN_SAT ⇒ **ε_sat = 1310 / 15 191 ≈ 0.0862**. Above ε_sat, increasing ε_org does **not** increase c_AN,bulk — it is fixed at 1310 mol/m³. The effect of higher ε_org enters entirely through D_i,mix (§4): organic droplets provide a parallel transport pathway for AN.

```julia
const M_AN             = 0.05306       # kg/mol
const RHO_AN           = 806.0         # kg/m³
const m_AN             = 11.59         # partition coefficient
const MOLAR_DENSITY_AN = RHO_AN / M_AN # ≈ 15 191 mol/m³
const C_AN_SAT         = MOLAR_DENSITY_AN / m_AN       # ≈ 1310 mol/m³
const EPS_ORG_SAT      = C_AN_SAT / MOLAR_DENSITY_AN   # ≈ 0.0862

function c_AN_bulk(eps_org)
    eps_org < EPS_ORG_SAT ? eps_org * MOLAR_DENSITY_AN : C_AN_SAT
end
```

> **ε_org = 0 is non-physical.** It represents zero AN in the system — no ADPN, no PN, no kinetic activity. All sweeps must start at ε_org > 0. The Stage 1 reference point is `ε_org = 0.02` (c_AN,bulk ≈ 304 mol/m³; single-phase).

Numerical table at the planned sweep points:

| ε_org | c_AN,bulk [mol/m³] | D_AN,mix [m²/s] | Regime |
|---|---|---|---|
| 0.02 | 304 | 2.30 × 10⁻⁹ | single-phase |
| 0.05 | 760 | 2.30 × 10⁻⁹ | single-phase |
| 0.08 | 1 215 | 2.30 × 10⁻⁹ | single-phase |
| **0.0862** | **1 310** (= C_AN_SAT) | 2.62 × 10⁻⁹ | threshold |
| 0.09 | 1 310 | 2.63 × 10⁻⁹ | two-phase |
| 0.15 | 1 310 | 2.86 × 10⁻⁹ | two-phase |
| 0.25 | 1 310 | 3.23 × 10⁻⁹ | two-phase |
| 0.30 | 1 310 | 3.41 × 10⁻⁹ | two-phase |

---

## 8. Three-Parameter Sweep

| Parameter | Range | Points |
|-----------|-------|--------|
| V (cathode potential) | −1.0 to −2.5 V vs SHE | Newton continuation |
| δ (boundary layer thickness) | 10, 20, 50, 100, 200 μm | 5 |
| ε_org (AN loading) | **0.02, 0.05, 0.08, 0.15, 0.25, 0.30** | 6 |

**The ε_org sweep spans both regimes:** {0.02, 0.05, 0.08} single-phase (below ε_sat ≈ 0.0862), {0.15, 0.25, 0.30} two-phase. The transition at ε_sat is directly observable. ε_org = 0 is never swept (non-physical).

**Total: 30 Newton continuation sweeps** → 3D performance map FE_ADPN(V, δ, ε_org).

The diffusion-only limiting current is:

$$j_{\mathrm{lim}} = \frac{n_e\,F\,D_{\mathrm{AN,mix}}\,c_{\mathrm{AN,bulk}}}{\delta}$$

Under the local-equilibrium assumption with the arithmetic-mean D_mix (two-phase regime), the enhancement at higher ε_org comes through the D_AN,mix > D_AN,aq channel. With the m_i-corrected D_eff (§4.2), the enhancement would be much larger. Comparing model predictions under both assumptions against Bloomquist data will determine which is correct.

---

## 9. Parameter Tables

### 9.1 Physical Constants

| Symbol | Value | Unit |
|--------|-------|------|
| F | 96,485.332 | C mol⁻¹ |
| R | 8.314463 | J mol⁻¹ K⁻¹ |
| T | 298.15 | K (25 °C) |

### 9.2 Kinetic Fitting Parameters

| Parameter | Initial | Range | Unit | Basis |
|-----------|---------|-------|------|-------|
| E⁰₁ (ADPN) | −1.3 | −1.0 to −1.5 | V vs SHE | Onset (Mathison JACS 2025) |
| E⁰₂ (PN) | −1.3 | −1.0 to −1.5 | V vs SHE | Same reaction centre |
| E⁰₃ (HER) | −0.83 | Fixed | V vs SHE | Nernst at pH 14 |
| j₀,₁ | 10⁻⁴ | 10⁻⁶ to 10⁻² | mA cm⁻² | Irreversible organic reduction |
| j₀,₂ | 10⁻⁴ | 10⁻⁶ to 10⁻² | mA cm⁻² | Similar to R1 |
| j₀,₃ (Cd) | 10⁻⁶ | 10⁻⁸ to 10⁻⁴ | mA cm⁻² | High HER overpotential |
| α_c,1 | 0.5 | 0.3–0.7 | — | Single-electron RDS |
| α_c,2 | 0.5 | 0.3–0.7 | — | Single-electron RDS |
| α_c,3 | 0.4 | 0.3–0.5 | — | Typical for Cd |

Multiply mA cm⁻² by 10 to convert to A m⁻² for use in the rate expressions.

### 9.3 Electrolyte Composition (Modestino Group)

| Component | Concentration | Role |
|-----------|--------------|------|
| Na₃PO₄ | 0.5 M (500 mol m⁻³) | Buffer + supporting electrolyte |
| TBA-OH | 0.02 M (20 mol m⁻³) | Selectivity agent — lumped into Na⁺ (see §2) |
| EDTA disodium | 0.03 M | Metal chelator — not tracked (negligible vs. phosphate) |
| Acrylonitrile | **ε_org ∈ [0.02, 0.30] volume fraction** | Reactant (sweep parameter; never 0) |
| Cathode | Cd foil | High hydrogen overpotential |

### 9.4 Partition Coefficients (for future m_i upgrade)

| Species | m_i = c_org/c_aq | Source |
|---------|------------------|--------|
| AN | 11.59 | Suwanvaipattana 2017 |
| SN | 7.72 | Suwanvaipattana 2017 |
| Trimer | 30 | Suwanvaipattana 2017 |
| Ether | 19.02 | Suwanvaipattana 2017 |
| All ions | 0 | Insoluble |

---

## 10. Numerical Methods

### 10.1 Mesh

Cell-centred finite-volume on a geometrically graded 1D mesh with **N_mesh = 100** cells. Cells are finest at x = 0 (electrode) and coarsen toward x = δ (bulk). Parameterised by stretch factor s = dx_max/dx_min (default s = 10):

```
r = s^(1/(N−1))              (common ratio)
dx_min = δ × (r−1) / (r^N − 1)
dx[k] = dx_min × r^(k−1)     for k = 1, …, N
```

At δ = 50 μm with N = 100 and s = 10: dx_min ≈ 0.13 μm, dx_max ≈ 1.3 μm. SG flux for charged species; centred differences for neutrals (reduces to SG with z = 0).

### 10.2 Jacobian — Block-Tridiagonal Bandwidth is **b = 17** (not 9)

The Jacobian has **block-tridiagonal structure** with block size 9 (one block per cell's DOFs). In cell-major layout, coupling of cell ix to cells ix−1, ix, ix+1 gives:

$$\text{half-bandwidth } b = 2 \cdot B - 1 = 2 \cdot 9 - 1 = \mathbf{17}$$

⚠️ **Correcting guide v4:** Earlier versions stated `b = 9`; that's the *block* size, not the bandwidth. With `b = 9` the banded-FD column-grouping at stride 19 aliases unrelated columns onto the same rows and produces an exactly-singular Jacobian. Use `b = 17` and `n_colors = 2b + 1 = 35` FD perturbations.

```julia
const JAC_BLOCK  = 9
const JAC_HALFBW = 2 * JAC_BLOCK - 1   # = 17
```

Banded FD: 35 residual evaluations per Jacobian (vs 900 for dense) → ~26× speedup per step. Sparse LU via `SparseArrays.jl`.

### 10.3 ForwardDiff AD Jacobian (Optional, Recommended for Hard Regions)

For regions where banded-FD truncation errors corrupt the Jacobian (cond(J) > 10¹⁵, high AN depletion, etc.), use `ForwardDiff.jacobian!` instead of FD. Requires:

1. **Type-generic residual chain.** All functions in the residual call graph must accept `AbstractVector{<:Real}` or `::Real` scalars, not `::Float64`. Work arrays must use `zeros(eltype(u), ...)` not `zeros(...)`.
2. **Taylor-smoothed Bernoulli branch** in `sg_flux` (§3.3). A hard `if abs(α) < 1e-10` branch gives a spurious zero derivative at the α=0 seed; the Taylor form `B(α) = 1 − α/2 + α²/12 − ...` for |α| < 0.01 is smooth everywhere.
3. **ForwardDiff.jl + ForwardDiff.JacobianConfig** with a chunk size (default 12).

```julia
ad_cfg = ForwardDiff.JacobianConfig(residual!, zeros(n), u, ForwardDiff.Chunk{12}())
J = zeros(n, n)
ForwardDiff.jacobian!(J, residual!, zeros(n), u, ad_cfg)
```

AD cost is ~n/chunk_size forward passes at ~2× Float64 cost each — about 4× slower per Jacobian than banded FD, but **exact to machine precision** and identical to FD in easy regions (verify with a comparison test at startup).

### 10.4 Newton Solver — Direct `(J + λI) du = −F` with Strict L2 Descent

Direct Newton with small Tikhonov damping `λ = 10⁻¹⁰` (keeps near-singular Jacobians solvable without distorting the Newton direction). Strict L2 descent `‖F_new‖₂ < ‖F_old‖₂` with 10-step backtracking line search. Step clamps per DOF type.

```julia
function newton_solve!(u, residual!;
                       tol = 1.0e-5, max_iter = 40,
                       lambda_fix = 1.0e-10,
                       jacobian_mode = :fd,              # :fd or :ad
                       max_backtrack = 10, verbose = false)
    # ... builds Jacobian (FD or AD), adds λ_fix to diagonal
    # ... solves (J + λI) du = -F
    # ... clamps du: log-conc ≤ 5.0, φ_ℓ ≤ 0.015 V
    # ... backtracks α starting from 1 until merit decreases strictly
end
```

| DOF type | Max step size per iteration |
|---|---|
| Log-concentration u_i | 5.0 (natural log units) |
| Electrolyte potential φ_l | 0.015 V |

Convergence tolerance: **‖F‖_∞ < 10⁻⁵**. This is loosened from 10⁻⁸ — charge conservation `|Σ z_i N_i|_interior < 10⁻¹⁰` is the real correctness gate (physicality check §13), and tight-tol wastes iterations chasing unreachable precision in near-singular regions.

> **Why not classical LM `(JᵀJ + λI)du = −JᵀF`?** Tested; heavy initial damping (`λ₀ = 1.0`) biased the continuation toward a non-physical solution branch. Direct `(J + λI)du = −F` with `λ = 10⁻¹⁰` preserves the Newton direction and stays on the physical attractor.

### 10.5 Continuation Strategy — V + Optional log-j Fallback

**Simple Newton V-continuation** (default). For each (δ, ε_org) pair, sweep V from −1.0 to −2.5 V vs SHE in steps of ds (initial 0.05 V), warm-starting Newton from the previous converged solution. Adaptive:

- Converges in ≤ 4 Newton iter → grow ds by 1.4× (up to ds_max = 0.20 V)
- Converges in 5–10 iter → grow by 1.1×
- Converges in 11–20 iter → keep ds
- Converges in > 20 iter → shrink by 0.7× (floor ds_min = 10⁻⁴ V)
- Fails (divergence, ‖F‖ not descending) → shrink by 0.3×; retry

**Optional log-j fallback** via `newton_continuation_logj`. When V-continuation stalls at a near-fold, switch to parameterisation by log₁₀(j_total) via augmented Newton with V as an extra DOF:

$$\text{augmented residual} = \begin{bmatrix} F(u, V) \\ \log_{10}(j_{\mathrm{tot}}(u, V)) - \log_{10} j_{\mathrm{target}} \end{bmatrix} = 0$$

Steps uniformly in log j. Equivalent to pseudo-arclength continuation with j as the arclength parameter. Primarily useful for dense sampling in the exponential-Tafel regime, not for extending the V-range (physical limits set the V floor).

### 10.6 Continuation Parameters α_buf and α_kin

Two scalar continuation parameters scale the source/flux terms during bootstrap (§12). Both are 0 at the start (residual is identically zero given the equilibrium initial guess) and ramped to 1 to recover the full physical model:

| Parameter | Scales | Effect at α = 0 | Effect at α = 1 |
|-----------|--------|-----------------|-----------------|
| α_buf | All buffer source terms R_buf,i | Buffer chemistry off — equilibrium IC is exact | Full buffer kinetics on |
| α_kin | All Faradaic flux densities j₁, j₂, j₃ | No current; flat profiles preserved | Full Tafel kinetics on |

### 10.7 Residual Form

> **Always use the integrated FV form** `F[ix] = J_left − J_right + S × dx` (residual in mol m⁻² s⁻¹), not the divided form `(J_right − J_left)/dx + S`. On the graded mesh (dx ratio up to 10:1 at s = 10), the divided form makes residual entries O(1/dx_min), degrading Jacobian row scaling.

### 10.8 Full Residual Assembly (Type-Generic)

DOF layout: cell-major with 9 DOFs per cell. Helper inlines:

```julia
@inline conc_dof(ix, k) = 9*(ix-1) + k    # k = 1..8 (log-concentration DOFs)
@inline phi_dof(ix)     = 9*ix             # one per cell
```

The residual signature is **type-generic** for ForwardDiff compatibility:

```julia
function full_residual!(res::AbstractVector{T}, u::AbstractVector{T},
                        mesh, eps_org::Float64, V::Real,
                        alpha_buf::Real, alpha_kin::Real,
                        c_eq) where {T<:Real}
    N = length(mesh.dx)
    c   = zeros(T, 8, N)
    phi = zeros(T, N)
    for ix in 1:N, k in 1:8
        c[k, ix] = exp(clamp(u[conc_dof(ix, k)], -50.0, 50.0))
    end
    for ix in 1:N
        phi[ix] = u[phi_dof(ix)]
    end

    # SG fluxes at interior faces (2..N)
    J = zeros(T, 8, N+1)
    for ix in 1:N-1
        dx_face = 0.5 * (mesh.dx[ix] + mesh.dx[ix+1])
        for k in 1:8
            J[k, ix+1] = sg_flux(c[k, ix], c[k, ix+1],
                                  phi[ix], phi[ix+1],
                                  D_mix(k, eps_org), z_species[k], dx_face)
        end
    end

    # Faradaic BC at face 1 (electrode, x = 0)
    j1, j2, j3 = tafel_currents(c[6, 1], phi[1], V, alpha_kin)
    J[2, 1] = +(j1 + j2 + j3) / F         # OH⁻
    J[6, 1] = -(2*j1 + j2) / (2*F)        # AN
    J[7, 1] = +(j1) / (2*F)               # ADPN
    J[8, 1] = +(j2) / (2*F)               # PN
    # J[1, 1], J[3:5, 1] = 0 (H⁺, phosphates: no Faradaic flux)

    R_buf = zeros(T, 8)
    for ix in 1:N
        if ix == N
            # Bulk Dirichlet (gauge: φ_l(δ) = 0)
            for k in 1:8
                res[conc_dof(ix, k)] = u[conc_dof(ix, k)] -
                                       log(bulk_concentration(k, c_eq, eps_org))
            end
            res[phi_dof(ix)] = phi[ix]
        else
            buffer_sources!(R_buf, c[1,ix], c[2,ix], c[3,ix], c[4,ix], c[5,ix], alpha_buf)
            dx_ix = mesh.dx[ix]
            for k in 1:8
                res[conc_dof(ix, k)] = J[k, ix] - J[k, ix+1] + R_buf[k] * dx_ix
            end
            # Current conservation: Σ z_k (J_left − J_right) over k = 1..5 (charged)
            s_L = 0.0; s_R = 0.0
            for k in 1:5
                s_L += z_species[k] * J[k, ix]
                s_R += z_species[k] * J[k, ix+1]
            end
            res[phi_dof(ix)] = s_L - s_R
        end
    end
    return res
end
```

Note the bulk AN concentration comes from `bulk_concentration(6, c_eq, eps_org)` which calls `c_AN_bulk(eps_org)` per §7.2.

---

## 11. Solution Caching

Binary files encoding ε_org, δ, and V in the filename: `s_eo0.150_d50_V-1.300000.bin`. Format: Int64 DOF count + Float64 vector. Cache immediately after Newton converges.

---

## 12. Implementation Stages (All Mandatory)

**Stage 1: NP + migration + buffer + kinetics at ε_org = 0.02 (single-phase reference).**

Minimal physical AN loading. D_i,mix = D_i,aq for all species (single-phase regime). Validates migration, SG scheme, buffer chemistry, and electrode kinetics. **Not ε_org = 0** — that's the pathological no-AN limit.

Sub-steps:
1. Compute `c_eq = solve_phosphate_equilibrium()` and build `u0 = make_initial_guess(N_mesh, c_eq, 0.02)`. At α_buf = 0, α_kin = 0, residual is identically zero.
2. **Buffer ramp:** α_buf from 0 → 1 in 10 uniform steps at V = −1.0 V vs SHE, α_kin = 0.
3. **Kinetics ramp:** α_kin geometric ×2 from 10⁻⁶ to 1.0 (21 steps). Still at V = −1.0 V.
4. **Newton continuation sweep over V** from −1.0 to −2.5 V.
5. Cache → plots → **STOP for review.**

**Stage 2: Activate D_i,mix(ε_org).**

Same equations, but D_i,mix now depends on ε_org across both regimes. Run at multiple ε_org values {0.02, 0.05, 0.08, 0.15, 0.25, 0.30}. The 0.08 → 0.15 transition directly probes the single-phase → two-phase transition. Compare polarisation curves and FE_ADPN to Stage 1 to quantify the D_mix effect.

Sub-steps: for each ε_org, warm-start from the Stage 1 converged solution at V = −1.0 V → Newton continuation sweep → cache → plots → **STOP.**

**Stage 3: Full 3D sweep.**

Runs Stage 2 across the complete (ε_org, δ) grid. 30 Newton continuation sweeps. Generate performance map and all comparison plots → **STOP.**

---

## 13. Physicality Checks

| Check | Expected |
|-------|----------|
| φ_l(x) profile | A few mV variation; larger at high j |
| Σ z_i J_i at each face | = 0 to machine precision |
| 9-species electroneutrality \|Σ z_i c_i\| | < 10⁻⁸ mol m⁻³ (trivially 0 by construction) |
| c_AN(x) at high j | Depletes toward electrode |
| Buffer residuals at bulk | < 10⁻⁶ mol m⁻³ s⁻¹ (Dirichlet exact) |
| pH(x) | Rises monotonically from bulk to surface at high j (OH⁻-pathway) |
| FE_ADPN vs ε_org | Peaks around ε_sat ≈ 0.086 ↑ in two-phase regime |
| D_AN,mix vs ε_org | Flat (single-phase) then rises (two-phase) |
| D_OH,mix vs ε_org | Flat (single-phase) then drops (two-phase) |
| No R_PT residual | Confirm R_PT is absent from residual |
| Bulk pH at x = δ | ≈ 13.03 (matches solve_phosphate_equilibrium with updated Ka) |

---

## 14. Module Structure

```
an_ehd/
├── params.jl           # Constants incl. MOLAR_DENSITY_AN, C_AN_SAT, EPS_ORG_SAT;
│                       # OH-pathway rate constants
├── mesh.jl             # make_mesh(N, delta; stretch)
├── diffusivity.jl      # D_mix(i, eps_org) — regime-aware (single vs two-phase)
├── chemistry.jl        # solve_phosphate_equilibrium, buffer_sources! (OH-pathway),
│                       # c_AN_bulk (Convention A), make_initial_guess
├── kinetics.jl         # j_ADPN, j_PN, j_HER (Tafel)
├── transport.jl        # sg_flux with Taylor-smoothed Bernoulli for |α| < 0.01
├── assembly.jl         # full_residual! (type-generic AbstractVector{T<:Real})
├── solver.jl           # newton_solve! (direct (J+λI)du=-F, :fd or :ad Jacobian);
│                       # newton_continuation; newton_continuation_logj
├── run_stage1.jl       # Stage 1 at ε_org = 0.02 → STOP
├── run_stage2.jl       # ε_org sweep {0.02, 0.05, 0.08, 0.15, 0.25, 0.30} → STOP
├── run_stage3.jl       # Full 3D sweep → STOP
├── plot_results.py     # 2×3 profile plots (incl. pH panel); 3-panel polarization
└── output/cache/
```

Note: `solve_phosphate_equilibrium` uses an inline bisection (not `Roots.jl`) to avoid Windows Defender Application Control blocking `Roots`'s DLL cache on some systems.

---

## 15. Common Pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| ε_org = 0 in sweep | c_AN,bulk = 0, log(0) = −Inf, NaN residual | Sweep from ε_org ≥ 0.02 |
| Using D_aq everywhere | ε_org has no effect past ε_sat | Use regime-aware D_mix |
| D_org ≠ 0 for ions | Unphysical ionic transport | Set D_org = 0 for all ions |
| Single-phase with mixed D | Spurious organic pathway when no droplets exist | D_mix = D_aq for ε_org < ε_sat |
| Adding ε_aq to R_buf | Double-counting | No prefactor (c_i is aqueous) |
| Including R_PT | Inconsistent with local equilibrium | Remove R_PT entirely |
| Wrong AN order | FE insensitive to j | c² for ADPN, c¹ for PN |
| H⁺-pathway buffer at high pH | Unphysical k_r ~ 10¹⁰ m³/(mol·s) | OH⁻-pathway stoichiometry (§6) |
| Cold start (no equilibrium IC) | Newton diverges immediately | solve_phosphate_equilibrium + make_initial_guess |
| SG overflow | exp crash | Clamp α to [−700, 700] |
| SG hard branch at α=0 | AD derivative spuriously 0 | Taylor expansion for \|α\| < 0.01 |
| Arithmetic D too weak | FE enhancement too small | Upgrade to m_i-corrected D_eff (§4.2) |
| Divided FV form on graded mesh | Jacobian ill-conditioned | Use integrated form J_L − J_R + S·dx (§10.7) |
| Missing φ_l gauge fix | Jacobian rank-deficient, Newton stalls | Set φ_l(δ) = 0 (Dirichlet, §7.2) |
| j₀ in mA/cm² passed directly | Rates off by factor of 10 | Convert: A/m² = mA/cm² × 10 |
| Dense FD Jacobian | ~60 ms/step, full sweep >10 min | Banded FD with 35-color grouping + sparse LU (§10.2) |
| Jacobian halfbandwidth b = 9 | Banded FD aliases columns, exactly singular | Use b = 2·block_size − 1 = 17 (§10.2) |
| Fixed absolute dx_min/dx_max in mesh | Domain length wrong across δ sweep | Use stretch-factor parameterisation |
| φ_l step clamp at 0.1 V | Over-large potential steps destabilise c_i DOFs | Clamp φ_l steps to 0.015 V |
| Implementing PAC before testing Newton cont. | Unnecessary complexity | Use simple Newton continuation first (§10.5) |
| Classical LM `(JᵀJ+λI)du=−JᵀF` with λ₀=1 | Picks wrong solution branch | Use direct `(J+λI)du=−F` with λ=10⁻¹⁰ (§10.4) |
| Merit slack 1.10 at tight tol | Accepts numerical-noise as progress | Strict L2 descent at tol ≤ 10⁻⁵ |
| Non-type-generic residual | ForwardDiff fails with type error | `AbstractVector{T<:Real}`; `zeros(eltype(u),...)` |

---

## 16. Potential Referencing

Model works internally in V vs SHE. Onsets from Mathison et al. (JACS 2025) on Cd: AN reduction −1.28 V; optimal ADPN −1.62 V; HER (no AN) −1.43 V.

| Reference electrode | Conversion |
|---------------------|------------|
| Ag/AgCl (sat. KCl) | E_SHE = E + 0.197 V |
| SCE | E_SHE = E + 0.241 V |
| RHE | E_SHE = E − 0.059 × pH (25°C) |

---

## 17. Experimental Data

**Bloomquist et al.** (CEJ 2026, 528, 172125): j = 70–300 mA cm⁻², ε_org = 0.05–0.30. FE_ADPN = 73–76% at j > 200 mA cm⁻². Bubble convection dominates. Note: experimentally-reported V span is 2–3× wider than the model's internal V-range because of uncompensated iR drop and TBA kinetic blocking (see §18).

**Mathison et al.** (JACS 2025, 147, 4296): Mechanism — radical coupling for ADPN, proton transfer for PN (strong KIE).

**Suwanvaipattana et al.** (J. Cleaner Prod. 2017, 142, 1296): D, m, d_p values, AN density (806 kg/m³).

**Huang et al.** (CEJ 2020, 382, 123006): Langmuir-adsorption kinetics on Pb.

**Costentin & Savéant** (J. Electroanal. Chem. 564, 2004, 99): ΔG⁰ = −1.84 eV for protonated radical coupling.

**Eigen & De Maeyer** (Z. Elektrochem. 1955): water autoprotolysis rate k₁,f = 1.4 × 10⁻³ M/s = 1.4 mol/(m³·s).

---

## 18. Fitting Strategy

1. **ε_org = 0.02, δ = 50 μm.** Fit kinetic parameters (j₀, α_c) to single-phase FE data.
2. **Activate D_mix at ε_org = 0.15.** Check whether FE_ADPN increases relative to single-phase. If not: the arithmetic D_mix is too weak → upgrade to m_i-corrected D_eff.
3. **δ sweep at ε_org = 0.15.** Match j_total trends to Bloomquist data.
4. **ε_org sweep at fixed δ.** Match FE vs organic loading curve.
5. **Full 3D sweep.** Generate performance map.

| Metric | Target |
|--------|--------|
| FE_ADPN peak (ε_org ≤ ε_sat, single-phase) | 50–60% |
| FE_ADPN peak (ε_org = 0.15, two-phase) | 73–80% |
| FE_ADPN at j = 200 mA cm⁻² | > 70% |
| FE enhancement from ε_org = 0.05 → 0.15 | +20–30 pp |

If the arithmetic D_mix cannot produce this enhancement, the m_i correction (§4.2) is needed. Optional physics additions — **kinetic saturation** `j_r = j_Tafel/(1 + j_Tafel/j_r,lim)` and **series resistance** `V_applied = V_interface + i·R_s` — allow direct V-axis overlay with Bloomquist cell voltages in a future extension, with j_r,lim and R_s as additional fit parameters.

---

## 19. Required Plots

Per-voltage profile figures use a **2×3 grid** to accommodate the pH panel:

| Panel | Content | Axes |
|---|---|---|
| (0,0) | H⁺ and OH⁻ (log y) | x [μm] / c [mol m⁻³] |
| (0,1) | Phosphate speciation | x [μm] / c [mol m⁻³] |
| (0,2) | **pH = −log₁₀(c_H/1000)** (with bulk reference line) | x [μm] / pH |
| (1,0) | AN / ADPN / PN | x [μm] / c [mol m⁻³] |
| (1,1) | φ_ℓ | x [μm] / mV |
| (1,2) | (reserved) | |

Per-sweep diagnostic panels:

| Panel | Content | Axes |
|---|---|---|
| (a1) | j_ADPN, j_PN, j_HER vs V — **log y** | −V vs SHE / mA cm⁻² |
| (a2) | j_ADPN, j_PN, j_HER vs V — **linear y** | −V vs SHE / mA cm⁻² |
| (b) | FE_ADPN, FE_PN, FE_HER vs V | −V vs SHE / % |
| (c) | FE_ADPN vs j at multiple ε_org | j / FE — key Bloomquist comparison |
| (d) | FE_ADPN vs ε_org at fixed j | ε_org / FE |
| (e) | c_AN(0)/c_AN,bulk vs j | j / depletion ratio |
| (f) | φ_l(0) vs j at multiple ε_org | j / mV (ohmic penalty) |
| (g) | D_AN,mix and D_OH,mix vs ε_org — **show regime transition at ε_sat** | ε_org / D |
| (h) | Production rate vs ε_org | ε_org / kg cm⁻² h⁻¹ |

Panel (g) is the most important diagnostic of the regime transition — it shows the step at ε_sat where single-phase D_aq gives way to the two-phase arithmetic mean.

---

*References: Bloomquist et al. CEJ 2026; Corpus et al. Joule 2023; Weng, Bell & Weber PCCP 2018; Huang et al. CEJ 2020; Suwanvaipattana et al. J. Cleaner Prod. 2017; Mathison et al. JACS 2025; Costentin & Savéant J. Electroanal. Chem. 2004; Lasia J. Electroanal. Chem. 1995; Eigen & De Maeyer Z. Elektrochem. 1955.*

*Guide v5 written 2026-04-21. Primary changes from v4: corrected AN bulk formula (Convention A, per total volume), regime-aware D_mix, OH⁻-pathway phosphate buffer chemistry with Eigen–De Maeyer water rates, ForwardDiff AD Jacobian option with Taylor-smoothed Bernoulli branch, direct `(J+λI)du=−F` Newton with tol=10⁻⁵, ε_org sweep shifted to [0.02, 0.30] (never 0). See CHANGELOG_V4toV5 for details.*
