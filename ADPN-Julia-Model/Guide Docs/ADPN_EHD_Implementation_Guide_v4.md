# Acrylonitrile Electrohydrodimerization (EHD) — 1D Planar Electrode Model

**Bui Lab | NYU Tandon School of Engineering | April 2026**

Scope: 1D Nernst diffusion layer (δ), Nernst–Planck transport with migration, Tafel kinetics for ADPN/PN/HER, phosphate buffer chemistry, local-equilibrium multiphase with arithmetic-mean D_i,mix, on a Cd cathode. **Operating temperature: T = 298.15 K (25 °C).**

---

## Convention

> All concentrations c_i in this guide are **aqueous-phase** values [mol m⁻³ of aqueous solution]. Under the local-equilibrium assumption (Da >> 1), organic and aqueous phases equilibrate instantaneously at every position: c_i,org(x) = m_i × c_i,aq(x). The organic phase acts as a parallel transport pathway captured entirely through the effective diffusivity D_i,mix. **No explicit phase transfer term R_PT appears** in the governing equations. No volume-fraction prefactors on any source terms.

---

## Table of Contents

1. Physical Domain
2. Species and Degrees of Freedom
3. Governing Equations
4. Mixture-Averaged Diffusivities
5. Electrochemical Kinetics
6. Phosphate Buffer Chemistry
   - 6.1 Equilibria
   - 6.2 R_buf,i Expressions
   - 6.3 Forward Rate Constants
   - 6.4 Bulk Equilibrium Initial Guess
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
 ADPN + PN + HER   |  D_i,mix(ε_org) eff. diffusivity   |  c_i = c_i,bulk
 Flux BCs           |  Phosphate buffer chemistry         |
                    |  No R_PT (local equilibrium)        |
                    |  Current conservation (φ_l)         |
```

Key findings from Bloomquist et al. (CEJ 2026, 528, 172125):

- FE_ADPN increases from <50% to >80% when ε_org exceeds the solubility limit (~0.09)
- Bubble-induced convection dominates over inlet flow regime
- High selectivity maintained at j > 200 mA cm⁻²
- ε_org explored from 0.05 to 0.30

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
| D_i,mix | Effective diffusivity through emulsion (§4) | m² s⁻¹ |
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

Standard centred-difference flux approximations become numerically unstable when the dimensionless migration Péclet number α = z_i F Δφ/(RT) exceeds ~2 per cell, producing spurious oscillations. The SG scheme uses an exponential fitting to handle arbitrarily large α stably. For neutral species (α = 0) it reduces to the standard centred difference:

$$J_i = -\frac{D_{i,\mathrm{mix}}}{\Delta x}\left[B(\alpha)\,c_R - B(-\alpha)\,c_L\right], \quad \alpha = \frac{z_i F(\phi_R - \phi_L)}{RT}$$

with the Bernoulli functions:

$$B(\alpha) = \frac{\alpha}{e^\alpha - 1}, \qquad B(-\alpha) = \frac{\alpha\,e^\alpha}{e^\alpha - 1}$$

where c_L, c_R are aqueous concentrations in the left and right cells [mol m⁻³]; φ_L, φ_R are electrolyte potentials [V]; and Δx is the distance between cell centres [m]. α is clamped to [−700, 700] to avoid floating-point overflow.

```julia
function sg_flux(c_L, c_R, phi_L, phi_R, D_mix, z_i, dx, T=298.15)
    z_i == 0 && return -D_mix * (c_R - c_L) / dx
    alpha = clamp(z_i * F * (phi_R - phi_L) / (R * T), -700.0, 700.0)
    abs(alpha) < 1e-10 && return -D_mix * (c_R - c_L) / dx
    B_pos = alpha / (exp(alpha) - 1.0)
    B_neg = alpha * exp(alpha) / (exp(alpha) - 1.0)
    return -D_mix / dx * (B_pos * c_R - B_neg * c_L)
end
```

---

## 4. Mixture-Averaged Diffusivities

### 4.1 Arithmetic Mean (Current Formulation)

The effective diffusivity is a simple volume-weighted arithmetic mean of the diffusivity in each phase:

$$D_{i,\mathrm{mix}} = \varepsilon_{\mathrm{org}}\,D_{i,\mathrm{org}} + (1 - \varepsilon_{\mathrm{org}})\,D_{i,\mathrm{aq}}$$

| Symbol | Definition | Unit |
|--------|-----------|------|
| ε_org | Organic-phase volume fraction (sweep parameter) | — |
| D_i,aq | Diffusivity of species i in pure aqueous phase | m² s⁻¹ |
| D_i,org | Diffusivity of species i in the organic phase | m² s⁻¹ |

This is the simplest mixing rule and does **not** include the partition coefficient m_i. The full local-equilibrium effective diffusivity is D_i,eff = (1−ε_org) D_i,aq + ε_org m_i D_i,org (see §4.2); the arithmetic mean is the first thing to test and will be upgraded if it underpredicts the experimental FE_ADPN enhancement.

Ions are insoluble in the organic phase, so D_i,org = 0 for all ions:

$$D_{i,\mathrm{mix}}^{\,\mathrm{ion}} = (1 - \varepsilon_{\mathrm{org}})\,D_{i,\mathrm{aq}}$$

This means increasing ε_org **reduces ionic transport**, creating an ohmic penalty. Neutral organic diffusivities are estimated via the Stokes–Einstein viscosity ratio: D_i,org ≈ D_i,aq × (μ_aq/μ_org) ≈ 2.6 × D_i,aq (μ_H₂O = 0.89 mPa·s, μ_AN = 0.34 mPa·s at 25°C).

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

D_mix(i, eps_org) = eps_org * D_org[i] + (1.0 - eps_org) * D_aq[i]
```

### Effect of ε_org on Mixture Diffusivities

| ε_org | D_AN,mix (m²/s) | D_OH,mix (m²/s) | D_AN/D_AN,aq | D_OH/D_OH,aq |
|-------|----------------|----------------|--------------|--------------|
| 0.00 | 2.30 × 10⁻⁹ | 5.27 × 10⁻⁹ | 1.00 | 1.00 |
| 0.09 | 2.63 × 10⁻⁹ | 4.80 × 10⁻⁹ | 1.14 | 0.91 |
| 0.15 | 2.86 × 10⁻⁹ | 4.48 × 10⁻⁹ | 1.24 | 0.85 |
| 0.25 | 3.23 × 10⁻⁹ | 3.95 × 10⁻⁹ | 1.40 | 0.75 |
| 0.30 | 3.41 × 10⁻⁹ | 3.69 × 10⁻⁹ | 1.48 | 0.70 |

AN transport improves by 48% at ε_org = 0.30 while OH⁻ transport degrades by 30% — this quantifies the selectivity–conductivity trade-off.

### 4.2 m_i-Corrected Upgrade (Future)

If the arithmetic mean is too weak to reproduce the experimental FE enhancement, the physically correct local-equilibrium effective diffusivity includes the partition coefficient m_i = c_i,org/c_i,aq:

$$D_{i,\mathrm{eff}} = (1-\varepsilon_{\mathrm{org}})\,D_{i,\mathrm{aq}} + \varepsilon_{\mathrm{org}}\,m_i\,D_{i,\mathrm{org}}$$

For AN (m_AN = 11.59) at ε_org = 0.25: D_AN,eff = 19.1 × 10⁻⁹ m²/s vs D_AN,mix = 3.23 × 10⁻⁹. For ions: m_i = 0, so D_ion,eff = (1−ε_org) D_i,aq — no change.

```julia
# One-line switch from arithmetic to m_i-corrected:
D_mix(i, eo) = eo * D_org[i] + (1-eo) * D_aq[i]          # arithmetic
D_eff(i, eo) = eo * m[i] * D_org[i] + (1-eo) * D_aq[i]   # m-corrected
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

```
j₁ = j₀,₁ × (c_AN / c_ref)² × exp(−α_c,1 × F × η₁ / (R×T))
```

**R2 — PN** (first-order in AN):

$$j_2 = j_{0,2}\,\frac{c_{\mathrm{AN}}}{c_{\mathrm{ref}}}\exp\!\left(-\frac{\alpha_{c,2} F\,\eta_2}{RT}\right)$$

```
j₂ = j₀,₂ × (c_AN / c_ref) × exp(−α_c,2 × F × η₂ / (R×T))
```

**R3 — HER** (independent of AN):

$$j_3 = j_{0,3}\exp\!\left(-\frac{\alpha_{c,3} F\,\eta_3}{RT}\right)$$

```
j₃ = j₀,₃ × exp(−α_c,3 × F × η₃ / (R×T))
```

The overpotential for each reaction:

$$\eta_r = (\phi_s - \phi_\ell) - E_r^{\,0}$$

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

> **Unit note:** j₀,r is given in the fitting table (§9.2) in mA cm⁻². Convert to SI via ×10: 1 mA cm⁻² = 10 A m⁻². c_ref = 1000 mol m⁻³ so that (c_AN/c_ref) is dimensionless; c_AN can exceed c_ref at high AN loading (above solubility, c_AN ≈ 1310 mol m⁻³), which is physical.

The c² vs c dependence is the central selectivity feature: as c_AN drops at the electrode under high current density, ADPN selectivity falls faster than PN.

The Faradaic efficiency for ADPN:

$$\mathrm{FE}_{\mathrm{ADPN}} = \frac{j_1}{j_1 + j_2 + j_3} \times 100\%$$

### 5.3 Alternative: Huang/Suwanvaipattana Rate Law

Huang et al. (CEJ 2020) use a Langmuir-adsorption formulation on Pb:

$$r_1 = k_1\left[\frac{K_{\mathrm{AN}}\,c_{\mathrm{AN}}}{1+K_{\mathrm{AN}}\,c_{\mathrm{AN}}}\right]^2, \qquad r_2 = k_2\,\frac{K_{\mathrm{AN}}\,c_{\mathrm{AN}}}{1+K_{\mathrm{AN}}\,c_{\mathrm{AN}}}$$

Fitted parameters (on Pb): k₁,ref = 0.56 ± 0.23 h⁻¹ m, k₂,ref = 0.075 ± 0.029 L mol⁻¹ h⁻¹ m, K_AN,ref = 10.80 ± 2.38 L mol⁻¹, ΔH_AN = −7.15 kJ mol⁻¹, α₁ = 1.48 V⁻¹, α₂ = 2.04 V⁻¹.

---

## 6. Phosphate Buffer Chemistry

### 6.1 Equilibria

Three equilibria are tracked. H₃PO₄ dissociation (pK = 2.12) is omitted because it is negligible at the operating pH of 12–13.

| # | Reaction | K (25°C, mol/m³ units) | pK |
|---|----------|------------------------|-----|
| 1 | H₂O ⇌ H⁺ + OH⁻ | K_w = 1.0 × 10⁻⁸ (mol/m³)² | 14.00 |
| 2 | H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻ | K_a2 = 6.2 × 10⁻⁵ mol/m³ | 7.21 |
| 3 | HPO₄²⁻ ⇌ H⁺ + PO₄³⁻ | K_a3 = 4.8 × 10⁻¹⁰ mol/m³ | 12.32 |

> **Unit convention:** All K values are in mol m⁻³ (SI). Relation to standard mol L⁻¹ values: K_w = 10⁻¹⁴ (mol/L)² × (1000)² = 10⁻⁸ (mol/m³)². The pK values listed are in standard units (−log₁₀[mol/L]).

Reverse rate constants are always derived: k_r = k_f / K_eq. Never set independently.

### 6.2 R_buf,i Expressions

For each equilibrium, the net rate [mol m⁻³ s⁻¹] is forward minus reverse:

$$r_1 = k_{1,f} - k_{1,r}\,c_{\mathrm{H}^+}\,c_{\mathrm{OH}^-}$$
$$r_2 = k_{2,f}\,c_{\mathrm{H_2PO_4^-}} - k_{2,r}\,c_{\mathrm{H}^+}\,c_{\mathrm{HPO_4^{2-}}}$$
$$r_3 = k_{3,f}\,c_{\mathrm{HPO_4^{2-}}} - k_{3,r}\,c_{\mathrm{H}^+}\,c_{\mathrm{PO_4^{3-}}}$$

Species allocation:

$$R_{\mathrm{buf}}[\mathrm{H}^+] = +r_1 + r_2 + r_3 \qquad R_{\mathrm{buf}}[\mathrm{OH}^-] = +r_1$$
$$R_{\mathrm{buf}}[\mathrm{H_2PO_4^-}] = -r_2 \qquad R_{\mathrm{buf}}[\mathrm{HPO_4^{2-}}] = +r_2 - r_3 \qquad R_{\mathrm{buf}}[\mathrm{PO_4^{3-}}] = +r_3$$

R_buf = 0 for AN, ADPN, and PN (no buffer reactions). **No prefactor** — c_i is aqueous; R_buf is rate per aqueous volume.

```julia
function buffer_sources!(R, c_H, c_OH, c_H2PO4, c_HPO4, c_PO4, alpha_buf=1.0)
    r1 = k1f - k1r * c_H * c_OH                   # H₂O ⇌ H⁺ + OH⁻
    r2 = k2f * c_H2PO4 - k2r * c_H * c_HPO4       # H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻
    r3 = k3f * c_HPO4  - k3r * c_H * c_PO4         # HPO₄²⁻ ⇌ H⁺ + PO₄³⁻
    R[1] = alpha_buf * (+r1 + r2 + r3)    # H⁺
    R[2] = alpha_buf * (+r1)              # OH⁻
    R[3] = alpha_buf * (-r2)              # H₂PO₄⁻
    R[4] = alpha_buf * (+r2 - r3)         # HPO₄²⁻
    R[5] = alpha_buf * (+r3)              # PO₄³⁻
    R[6:8] .= 0.0                          # AN, ADPN, PN
end
```

The optional `alpha_buf` argument (default 1.0) is the buffer continuation parameter used during bootstrap (see §10.6, §12).

### 6.3 Forward Rate Constants

| k | Value | Unit |
|---|-------|------|
| k₁,f | 2.4 × 10⁻⁵ | s⁻¹ |
| k₂,f | ~10⁶ | s⁻¹ |
| k₃,f | ~10² | s⁻¹ |

### 6.4 Bulk Equilibrium Initial Guess

**This is the most critical step for convergence.** Buffer reactions are stiff: forward and reverse rates are both large but nearly cancel at equilibrium. An initial guess that is even slightly off-equilibrium produces buffer source terms of O(10³–10⁵) mol m⁻³ s⁻¹ — orders of magnitude larger than transport fluxes — causing Newton to fail or stall for hundreds of iterations even at α_buf = 0.001.

**The fix: solve the phosphate charge balance exactly before starting Newton.**

The equilibrium speciation follows from two constraints:

**Phosphate conservation:**
$$c_{\mathrm{H_2PO_4}} + c_{\mathrm{HPO_4}} + c_{\mathrm{PO_4}} = C_{P,\mathrm{total}}$$

Expressing each species in terms of c_H alone:
$$c_{\mathrm{HPO_4}} = \frac{C_{P,\mathrm{total}}}{c_H/K_{a2} + 1 + K_{a3}/c_H}, \quad c_{\mathrm{H_2PO_4}} = \frac{c_H \, c_{\mathrm{HPO_4}}}{K_{a2}}, \quad c_{\mathrm{PO_4}} = \frac{K_{a3} \, c_{\mathrm{HPO_4}}}{c_H}$$

**Charge balance** (one equation in c_H):
$$c_{\mathrm{Na,input}} + c_H = \frac{K_w}{c_H} + c_{\mathrm{H_2PO_4}} + 2\,c_{\mathrm{HPO_4}} + 3\,c_{\mathrm{PO_4}}$$

where c_Na,input = 3 × C_Na₃PO₄ + C_TBA-OH = 3 × 500 + 20 = **1520 mol m⁻³** (TBA⁺ lumped into Na⁺, see §2).

This is a 1D root-finding problem in log₁₀(c_H):

```julia
const K_w  = 1.0e-8    # (mol/m³)², water autoprotolysis at 25°C
const K_a2 = 6.2e-5    # mol/m³, H₂PO₄⁻ ↔ H⁺ + HPO₄²⁻ at 25°C
const K_a3 = 4.8e-10   # mol/m³, HPO₄²⁻ ↔ H⁺ + PO₄³⁻ at 25°C

"""
    solve_phosphate_equilibrium(C_P_total, c_Na_input) → NamedTuple

Solve the phosphate buffer charge balance in the well-mixed bulk.
All concentrations in mol m⁻³.

  C_P_total  : total phosphate [mol/m³]  (0.5 M Na₃PO₄ → 500)
  c_Na_input : total monovalent cation input [mol/m³]
               (3×500 Na⁺ from Na₃PO₄ + 20 TBA⁺ from TBA-OH → 1520)

Returns NamedTuple (H, OH, H2PO4, HPO4, PO4, Na, pH).
"""
function solve_phosphate_equilibrium(C_P_total=500.0, c_Na_input=1520.0)

    function charge_balance(log10_cH)
        cH    = 10.0^log10_cH
        cOH   = K_w / cH
        denom = cH / K_a2 + 1.0 + K_a3 / cH
        cHPO4  = C_P_total / denom
        cH2PO4 = cH * cHPO4 / K_a2
        cPO4   = K_a3 * cHPO4 / cH
        # charge balance: (fixed cations + H⁺) − anions = 0
        return c_Na_input + cH - (cOH + cH2PO4 + 2.0*cHPO4 + 3.0*cPO4)
    end

    # Bisect over pH 6–15 in standard units
    # → log₁₀(c_H [mol/m³]) ∈ [−12, −3]
    using Roots
    log10_cH_eq = find_zero(charge_balance, (-12.0, -3.0), Bisection())

    cH     = 10.0^log10_cH_eq
    cOH    = K_w / cH
    denom  = cH / K_a2 + 1.0 + K_a3 / cH
    cHPO4  = C_P_total / denom
    cH2PO4 = cH * cHPO4 / K_a2
    cPO4   = K_a3 * cHPO4 / cH
    pH_std = -log10(cH / 1000.0)   # pH in standard mol L⁻¹ units

    return (H=cH, OH=cOH, H2PO4=cH2PO4, HPO4=cHPO4, PO4=cPO4,
            Na=c_Na_input, pH=pH_std)
end
```

**Verification — always run before starting Newton:**

```julia
c_eq = solve_phosphate_equilibrium()

# Equilibrium product checks
@assert abs(c_eq.H * c_eq.OH              - K_w)  / K_w  < 1e-8   "K_w mismatch"
@assert abs(c_eq.H * c_eq.HPO4 / c_eq.H2PO4 - K_a2) / K_a2 < 1e-8   "K_a2 mismatch"
@assert abs(c_eq.H * c_eq.PO4  / c_eq.HPO4  - K_a3) / K_a3 < 1e-8   "K_a3 mismatch"

# Buffer sources must be zero at equilibrium
R_check = zeros(8)
buffer_sources!(R_check, c_eq.H, c_eq.OH, c_eq.H2PO4, c_eq.HPO4, c_eq.PO4)
@assert maximum(abs.(R_check)) < 1e-6  "Buffer not at equilibrium: max|R_buf| = $(maximum(abs.(R_check))) mol/m³/s"
```

**Expected equilibrium values** for 0.5 M Na₃PO₄ + 0.02 M TBA-OH at 25°C:

| Species | Approximate value | Unit |
|---------|-------------------|------|
| c_H | ~10⁻¹⁰ | mol m⁻³ |
| c_OH | ~100 | mol m⁻³ |
| c_H₂PO₄ | ~0 | mol m⁻³ |
| c_HPO₄ | ~70–90 | mol m⁻³ |
| c_PO₄ | ~410–430 | mol m⁻³ |
| c_Na | 1520 | mol m⁻³ |
| pH | ~13.0 | — |

Most phosphate sits as PO₄³⁻ since pH ≈ 13.0 is well above pKa3 = 12.32.

**Building the full DOF initial guess vector:**

```julia
# DOF layout: 8 log-concentrations + 1 φ_l per cell, cell-major (interleaved)
# Species order: H⁺(1), OH⁻(2), H₂PO₄⁻(3), HPO₄²⁻(4), PO₄³⁻(5), AN(6), ADPN(7), PN(8)
function make_initial_guess(N_mesh, c_eq, c_AN_bulk; c_seed=1e-3)
    u0 = zeros(9 * N_mesh)
    log_c0 = [log(c_eq.H),    log(c_eq.OH),
              log(c_eq.H2PO4), log(c_eq.HPO4), log(c_eq.PO4),
              log(c_AN_bulk),
              log(c_seed),     # ADPN: small seed to avoid log(0)
              log(c_seed)]     # PN:   small seed to avoid log(0)
    for ix in 1:N_mesh
        u0[9*(ix-1)+1 : 9*(ix-1)+8] .= log_c0
        u0[9*(ix-1)+9] = 0.0   # φ_l = 0 everywhere
    end
    return u0
end
```

With this initialisation and α_buf = 0 (buffer off), α_kin = 0 (kinetics off), the residual is **identically zero**: no concentration gradients, no current, no buffer sources. Newton converges in 0 iterations — the ideal starting point for the bootstrap protocol (§12).

---

## 7. Boundary Conditions

### 7.1 Electrode (x = 0): Faradaic Flux

The species fluxes at the electrode surface are set by Faraday's law applied to the three cathodic reactions. Two AN molecules are consumed per ADPN produced and one AN per PN produced, giving:

$$N_{\mathrm{AN}}\big|_{x=0} = -\frac{2j_1 + j_2}{2F}$$

ADPN and PN are produced at the surface and diffuse away:

$$N_{\mathrm{ADPN}}\big|_{x=0} = +\frac{j_1}{2F}$$

$$N_{\mathrm{PN}}\big|_{x=0} = +\frac{j_2}{2F}$$

Hydroxide is produced by all three reactions (2 OH⁻ per 2 electrons transferred):

$$N_{\mathrm{OH}^-}\big|_{x=0} = +\frac{j_1 + j_2 + j_3}{F}$$

All other species are non-electroactive and have zero flux:

$$N_i\big|_{x=0} = 0 \quad \text{for } i \in \{\mathrm{H}^+,\ \text{phosphates}\}$$

The electrolyte potential φ_l(0) is **not prescribed** — it is determined by the current conservation equation (§3.2). The overpotentials in the Tafel rate expressions use this solved value.

### 7.2 Bulk (x = δ): Dirichlet

At the bulk boundary, all species are fixed at their well-mixed bulk values and the electrolyte potential is referenced to zero:

$$c_i(\delta) = c_{i,\mathrm{bulk}} \qquad \phi_\ell(\delta) = 0$$

The bulk equilibrium concentrations of H⁺, OH⁻, and the three phosphate species are computed by `solve_phosphate_equilibrium()` (§6.4). The φ_l(δ) = 0 Dirichlet also serves as the gauge fix required for well-posedness of the current conservation equation (see §3.2).

ADPN and PN bulk concentrations are set to small seed values (~10⁻³ mol m⁻³) to avoid log(0) in the initial guess — in an ideal well-mixed bulk they would be zero.

The bulk aqueous AN concentration c_AN,bulk depends on ε_org because the organic phase sequesters AN. The total amount of AN added to the system partitions between the aqueous phase (concentration c_AN,aq) and the organic phase (concentration c_AN,org = m_AN × c_AN,aq). Two regimes arise:

**Below solubility** (ε_org < 0.09): All AN dissolves in the aqueous phase. The aqueous volume is reduced by the organic fraction, so the dissolved concentration is the total AN divided by the aqueous volume fraction:

$$c_{\mathrm{AN,bulk}} = \frac{C_{\mathrm{AN,total}}}{1 - \varepsilon_{\mathrm{org}}}$$

where C_AN,total is the moles of AN per unit total volume.

**Above solubility** (ε_org ≥ 0.09): The aqueous phase is saturated. Additional AN remains in the organic phase as undissolved droplets. The aqueous concentration is pinned at the equilibrium solubility set by the partition coefficient:

$$c_{\mathrm{AN,bulk}} = c_{\mathrm{AN,eq}} = \frac{\rho_{\mathrm{AN}}}{M_{\mathrm{AN}}\,m_{\mathrm{AN}}} \approx 1{,}310\ \mathrm{mol\,m^{-3}}$$

Above the solubility limit, increasing ε_org does **not** increase c_AN,bulk — it is fixed at ~1,310 mol m⁻³. The effect of higher ε_org enters entirely through D_i,mix: the organic droplets provide a parallel transport pathway for AN, increasing the effective flux of AN toward the electrode without changing the boundary concentration.

---

## 8. Three-Parameter Sweep

| Parameter | Range | Points |
|-----------|-------|--------|
| V (cathode potential) | −1.0 to −2.5 V vs SHE | ~100 (Newton continuation) |
| δ (boundary layer thickness) | 10, 20, 50, 100, 200 μm | 5 |
| ε_org (organic volume fraction) | 0.00, 0.05, 0.09, 0.15, 0.25, 0.30 | 6 |

**Total: 30 Newton continuation sweeps** → 3D performance map FE_ADPN(V, δ, ε_org).

The diffusion-only limiting current is:

$$j_{\mathrm{lim}} = \frac{n_e\,F\,D_{\mathrm{AN,mix}}\,c_{\mathrm{AN,bulk}}}{\delta}$$

Under the local-equilibrium assumption with the arithmetic-mean D_mix, the enhancement comes purely from D_AN,mix > D_AN,aq. With the m_i-corrected D_eff (§4.2), the enhancement would be much larger. Comparing model predictions under both assumptions against Bloomquist data will determine which is correct.

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
| Acrylonitrile | 0.05–0.30 vol fraction | Reactant (sweep parameter) |
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

Cell-centred finite-volume on a geometrically graded 1D mesh with **N_mesh = 100** cells. Cells are finest at x = 0 (electrode) and coarsen toward x = δ (bulk). The mesh is parameterised by a **stretch factor** s = dx_max/dx_min (default s = 10). Given s and N, the minimum cell width is derived analytically from the geometric series sum:

```
r = s^(1/(N−1))              (common ratio)
dx_min = δ × (r−1) / (r^N − 1)
```

Cell widths: `dx[k] = dx_min × r^(k−1)` for k = 1, …, N. This construction guarantees the mesh fills [0, δ] exactly for any δ — a fixed absolute dx_min would not.

```julia
function make_mesh(N::Int, delta::Float64; stretch::Float64=10.0)
    r      = stretch^(1.0 / (N - 1))
    dx_min = delta * (r - 1.0) / (r^N - 1.0)
    dx     = [dx_min * r^(k-1) for k in 1:N]
    # Cell-centre positions
    x_c    = cumsum(dx) .- dx ./ 2.0
    return (dx=dx, x_c=x_c)
end
```

At δ = 50 μm with N = 100 and s = 10: dx_min ≈ 0.10 μm, dx_max ≈ 1.0 μm, dx ratio = 10. Across the full δ sweep (10–200 μm), dx_min scales proportionally — resolution is always δ/100 at coarsest and δ/1000 at finest. SG flux for charged species; centred differences for neutrals.

### 10.2 Jacobian Strategy

The Jacobian has **block-tridiagonal structure** with block size 9 (one per DOF) and half-bandwidth b = 9 (lower = upper). Total bandwidth = 2b = 18.

**Banded finite-difference with column coloring.** The banded structure means that columns separated by more than 2b+1 = 19 positions never share a non-zero row — they can be perturbed simultaneously. A sequential-column-grouping coloring therefore requires exactly **19 perturbation vectors** to reconstruct the full Jacobian, regardless of N:

```julia
# Banded FD Jacobian via column coloring
# bandwidth b = 9 (block size), n_colors = 2b+1 = 19
function banded_fd_jacobian!(J_sparse, residual!, u, F0, b=9)
    n = length(u)
    n_colors = 2*b + 1
    for color in 1:n_colors
        du = zeros(n)
        for j in color:n_colors:n
            du[j] = max(1e-7 * abs(u[j]), 1e-9)
        end
        dF = (residual!(u .+ du) .- F0) ./ du[color:n_colors:n]
        # scatter dF back into the appropriate diagonals of J_sparse
        # (store J_sparse as SparseMatrixCSC or BandedMatrix)
    end
end
```

Jacobian cost: 19 residual evaluations per Newton step vs 900 for dense FD → **47× speedup per step**. Solve cost: sparse/banded LU factorisation (O(N × b²)) vs dense LU (O(N³)) → further ~50× at N = 100. Use `SparseArrays.jl` with `LinearAlgebra.lu` or `KLU.jl` for factorisation.

> **Perturbation sizes:** Use ε_conc = 10⁻⁷ × max(|u_i|, 1) for log-concentration DOFs; ε_phi = 10⁻⁵ V for φ_l DOFs. Perturbing u_i = ln(c_i) corresponds to a multiplicative perturbation of c_i — this is more uniform across the wide concentration range (10⁻¹⁰ to 1520 mol m⁻³) than perturbing c_i directly.

### 10.3 Newton Solver

Newton with row scaling (divide each residual row by its L∞ norm at the first iteration) and step clamping:

| DOF type | Max step size |
|----------|--------------|
| Log-concentration u_i | 5.0 (natural log units) |
| Electrolyte potential φ_l | **0.015 V** |

The φ_l clamp is tightened to 0.015 V (from 0.1 V in v3) because the physical φ_l variation across the diffusion layer is only a few mV. A 0.1 V clamp allows Newton to take unphysically large potential steps that couple back into the SG fluxes and destabilise the concentration DOFs.

Convergence tolerance: ||F||_∞ < 10⁻⁴.

### 10.4 Continuation Strategy

**Simple Newton continuation** (no pseudo-arclength). For each (δ, ε_org) pair, sweep V from −1.0 to −2.5 V vs SHE in uniform steps of ds = 0.05 V, warm-starting Newton from the previous converged solution:

- If Newton converges in ≤ 5 iterations: increase ds by factor 1.5 (up to ds_max = 0.20 V).
- If Newton fails (diverges or exceeds 25 iterations): halve ds and retry from the last converged point.
- If ds < ds_min = 0.005 V: report convergence failure and skip this (V, δ, ε_org) point.

PAC (pseudo-arclength continuation) adds significant implementation complexity (augmented system, arclength parameterisation, tangent predictor) for a system with no expected fold bifurcations. This 1D NP diffusion layer is monotone in V — there is no turning point or bistability. Add PAC later only if the mass-transfer-limited regime shows repeated convergence failures with simple continuation.

### 10.5 Residual Form

> **Always use the integrated FV form** `F[ix] = J_left − J_right + S × dx` (residual in mol m⁻² s⁻¹), not the divided form `(J_right − J_left)/dx + S`. On the graded mesh (dx ratio up to 10:1 at s = 10), the divided form makes residual entries O(1/dx_min), degrading Jacobian row scaling.

### 10.6 Continuation Parameters α_buf and α_kin

Two scalar continuation parameters scale the source/flux terms during bootstrap (§12). Both are 0 at the start (residual is identically zero given the equilibrium initial guess) and ramped to 1 to recover the full physical model:

| Parameter | Scales | Effect at α = 0 | Effect at α = 1 |
|-----------|--------|-----------------|-----------------|
| α_buf | All buffer source terms R_buf,i | Buffer chemistry off — equilibrium IC is exact | Full buffer kinetics on |
| α_kin | All Faradaic flux densities j₁, j₂, j₃ | No current; flat profiles preserved | Full Tafel kinetics on |

Implementation: pass α_buf and α_kin into the buffer and kinetics evaluators as multiplicative scalars on the rate terms.

```julia
function buffer_sources!(R, c_H, c_OH, c_H2PO4, c_HPO4, c_PO4, alpha_buf=1.0)
    r1 = k1f - k1r * c_H * c_OH
    r2 = k2f * c_H2PO4 - k2r * c_H * c_HPO4
    r3 = k3f * c_HPO4  - k3r * c_H * c_PO4
    R[1] = alpha_buf * (+r1 + r2 + r3)
    R[2] = alpha_buf * (+r1)
    R[3] = alpha_buf * (-r2)
    R[4] = alpha_buf * (+r2 - r3)
    R[5] = alpha_buf * (+r3)
    R[6:8] .= 0.0
end

function tafel_currents(c_AN_surface, phi_l_surface, phi_s, alpha_kin=1.0)
    eta1 = (phi_s - phi_l_surface) - E0_1
    eta2 = (phi_s - phi_l_surface) - E0_2
    eta3 = (phi_s - phi_l_surface) - E0_3
    j1 = j0_1 * (c_AN_surface / c_ref)^2 * exp(-alpha_c1 * F * eta1 / (R * T))
    j2 = j0_2 * (c_AN_surface / c_ref)    * exp(-alpha_c2 * F * eta2 / (R * T))
    j3 = j0_3 *                              exp(-alpha_c3 * F * eta3 / (R * T))
    return alpha_kin * j1, alpha_kin * j2, alpha_kin * j3
end
```

### 10.7 Full Residual Assembly

DOF layout: cell-major with 9 DOFs per cell. For cell ix and species k ∈ {1..8}, the log-concentration index is `9*(ix−1) + k`; the φ_l index is `9*ix`. Helper inlines:

```julia
@inline conc_dof(ix, k) = 9*(ix-1) + k    # k = 1..8 (log-concentration DOFs)
@inline phi_dof(ix)     = 9*ix             # one per cell
```

Each cell contributes 9 residual rows: 8 species transport equations + 1 current conservation equation. Faces are indexed 1..N+1, where face 1 is the electrode (x = 0), face ix+1 is between cells ix and ix+1, and face N+1 is the bulk (x = δ). The first cell's left-face flux is the Faradaic boundary flux (§7.1); the last cell's residual is replaced by Dirichlet equations (§7.2).

```julia
function full_residual!(F, u, mesh, eps_org, V, alpha_buf, alpha_kin, c_eq)
    N = length(mesh.dx)
    z = (+1, -1, -1, -2, -3, 0, 0, 0)   # H⁺ OH⁻ H₂PO₄⁻ HPO₄²⁻ PO₄³⁻ AN ADPN PN

    # ---- Decode DOFs ----
    c   = zeros(8, N)
    phi = zeros(N)
    for ix in 1:N
        for k in 1:8
            c[k, ix] = exp(clamp(u[conc_dof(ix, k)], -50.0, 50.0))
        end
        phi[ix] = u[phi_dof(ix)]
    end

    # ---- Faces 2..N: SG fluxes between adjacent cells ----
    J = zeros(8, N+1)
    for ix in 1:N-1
        dx_face = 0.5*(mesh.dx[ix] + mesh.dx[ix+1])
        for k in 1:8
            D_k = D_mix(k, eps_org)
            J[k, ix+1] = sg_flux(c[k, ix], c[k, ix+1],
                                 phi[ix], phi[ix+1], D_k, z[k], dx_face)
        end
    end

    # ---- Face 1 (x = 0): Faradaic boundary fluxes ----
    # Sign: N_i is the +x flux; consumption at electrode → N_i(0) < 0,
    #       production at electrode → N_i(0) > 0.
    j1, j2, j3 = tafel_currents(c[6, 1], phi[1], V, alpha_kin)
    J[1, 1] = 0.0                            # H⁺   (no Faradaic flux)
    J[2, 1] = +(j1 + j2 + j3) / F            # OH⁻ produced
    J[3, 1] = 0.0; J[4, 1] = 0.0; J[5, 1] = 0.0   # phosphates: no flux
    J[6, 1] = -(2*j1 + j2) / (2*F)           # AN consumed
    J[7, 1] = +(j1) / (2*F)                  # ADPN produced
    J[8, 1] = +(j2) / (2*F)                  # PN produced
    # Face N+1 (x = δ) is unused — Dirichlet override below replaces last-cell residual.

    # ---- Assemble per-cell residuals ----
    R_buf = zeros(8)
    for ix in 1:N
        if ix == N
            # Dirichlet override at the bulk cell
            for k in 1:8
                F[conc_dof(ix, k)] = u[conc_dof(ix, k)] - log(bulk_concentration(k, c_eq, eps_org))
            end
            F[phi_dof(ix)] = phi[ix]              # gauge fix φ_l(δ) = 0
        else
            buffer_sources!(R_buf, c[1,ix], c[2,ix], c[3,ix], c[4,ix], c[5,ix], alpha_buf)
            # Species: integrated FV form  J_left − J_right + S·dx
            for k in 1:8
                F[conc_dof(ix, k)] = J[k, ix] - J[k, ix+1] + R_buf[k] * mesh.dx[ix]
            end
            # Current conservation: Σ z_i (J_left − J_right) over charged species (k = 1..5)
            sum_left  = 0.0
            sum_right = 0.0
            for k in 1:5
                sum_left  += z[k] * J[k, ix]
                sum_right += z[k] * J[k, ix+1]
            end
            F[phi_dof(ix)] = sum_left - sum_right
        end
    end
    return F
end
```

The Tafel rate evaluation uses `c[6, 1]` (cell-centred AN at the first cell) as the surface concentration approximation, and the *solved* `phi[1]` as the φ_l at the electrode for the overpotential `η = (V − φ_l) − E⁰`. This couples the Faradaic BC implicitly to the interior solution.

```julia
function bulk_concentration(k, c_eq, eps_org)
    k == 1 && return c_eq.H
    k == 2 && return c_eq.OH
    k == 3 && return c_eq.H2PO4
    k == 4 && return c_eq.HPO4
    k == 5 && return c_eq.PO4
    k == 6 && return c_AN_bulk(eps_org)
    k == 7 && return 1e-3                       # ADPN seed
    k == 8 && return 1e-3                       # PN seed
end

const C_AN_SAT = 1310.0        # mol m⁻³, aqueous saturation concentration of AN
function c_AN_bulk(eps_org; eps_sat=0.09)
    eps_org < eps_sat ? C_AN_SAT / (1.0 - eps_org) : C_AN_SAT
end
```

### 10.8 Newton Solver

```julia
using SparseArrays, LinearAlgebra

function newton_solve!(u, residual!; tol=1e-4, max_iter=25, verbose=false)
    n  = length(u)
    F0 = zeros(n)
    J  = build_banded_pattern(n, 9)        # block-tridiagonal sparsity, block size 9

    for iter in 1:max_iter
        residual!(F0, u)
        normF = maximum(abs.(F0))
        verbose && @info "Newton" iter normF
        normF < tol && return (converged=true, iter=iter, normF=normF)

        banded_fd_jacobian!(J, residual!, u, F0)   # 19 perturbations (§10.2)
        du = -(J \ F0)                              # sparse LU solve

        # Per-DOF-type step clamping (§10.3)
        N_cells = n ÷ 9
        for ix in 1:N_cells
            for k in 1:8
                idx = conc_dof(ix, k)
                du[idx] = clamp(du[idx], -5.0, 5.0)
            end
            idx_phi = phi_dof(ix)
            du[idx_phi] = clamp(du[idx_phi], -0.015, 0.015)
        end
        u .+= du
    end
    return (converged=false, iter=max_iter, normF=maximum(abs.(F0)))
end

function newton_continuation(u0, V_start, V_end, build_residual;
                              ds_init=0.05, ds_min=0.005, ds_max=0.20)
    u  = copy(u0)
    V  = V_start
    ds = ds_init
    history = Tuple{Float64,Vector{Float64}}[]
    while V ≥ V_end - 1e-12                   # V_end is more negative than V_start
        residual! = build_residual(V)         # closure capturing mesh, ε_org, δ, α_*, c_eq
        u_trial = copy(u)
        result = newton_solve!(u_trial, residual!)
        if result.converged
            push!(history, (V, copy(u_trial)))
            u  = u_trial
            ds = result.iter ≤ 5 ? min(ds * 1.5, ds_max) : ds
            V -= ds                            # step toward more negative V
        else
            ds /= 2.0
            ds < ds_min && error("Continuation failed at V = $V")
            # do not advance; retry with smaller ds
        end
    end
    return history
end
```

Bootstrap (Stage 1, §12) calls `newton_solve!` with progressively larger (α_buf, α_kin) values. Once both reach 1.0, `newton_continuation` is invoked to sweep V. Each continuation step warm-starts from the previous converged solution.

> **Banded sparsity pattern:** `build_banded_pattern(n, 9)` constructs a `SparseMatrixCSC` with non-zeros at every (i, j) where |i − j| ≤ 9 (the block-tridiagonal structure of block size 9). `banded_fd_jacobian!` then writes the 19 column-grouped FD perturbations into the existing pattern's non-zeros, avoiding repeated allocation.

---

## 11. Solution Caching

Binary files encoding ε_org, δ, and V in the filename: `s_eo0.150_d50_V-1.300000.bin`. Format: Int64 DOF count + Float64 vector. Cache immediately after Newton converges.

```julia
function save_solution(u::Vector{Float64}, eps_org, delta_um, V)
    fn = @sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, round(Int, delta_um*1e6), V)
    open(joinpath(CACHE_DIR, fn), "w") do f
        write(f, Int64(length(u)))
        write(f, u)
    end
end

function load_solution(eps_org, delta_um, V)
    fn = @sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, round(Int, delta_um*1e6), V)
    path = joinpath(CACHE_DIR, fn)
    isfile(path) || return nothing
    open(path, "r") do f
        n = read(f, Int64)
        read!(f, Vector{Float64}(undef, n))
    end
end
```

---

## 12. Implementation Stages (All Mandatory)

**Stage 1: NP + migration + buffer + kinetics at ε_org = 0.**

Single-phase reference case. No organic droplets, D_i,mix = D_i,aq. Validates migration, SG scheme, buffer chemistry, and electrode kinetics.

Sub-steps:
1. Compute `c_eq = solve_phosphate_equilibrium()` and build `u0 = make_initial_guess(...)` (§6.4). At α_buf = 0, α_kin = 0, residual is identically zero — Newton converges in 0 iterations.
2. **Buffer ramp:** increase α_buf from 0 → 1 in **10 uniform steps** (Δα = 0.1) at V = −1.0 V vs SHE, with α_kin = 0. Each step warm-starts from the previous converged solution. Buffer chemistry is smooth (Da >> 1, no thin reaction zones), so 10 steps is more than sufficient.
3. **Kinetics ramp:** increase α_kin geometrically — multiply by 2 each step — starting at α_kin = 10⁻⁶ and ending at α_kin = 1.0 (about **20 steps**: 10⁻⁶, 2×10⁻⁶, 4×10⁻⁶, …, 0.5, 1.0). Still at V = −1.0 V. The geometric ramp handles the large dynamic range of kinetic fluxes without requiring hundreds of steps.
4. **Newton continuation sweep over V** from −1.0 to −2.5 V (see §10.4 for step-size control).
5. Cache → plots → **STOP for review.**

**Stage 2: Activate D_i,mix(ε_org).**

Same equations, but D_i,mix now depends on ε_org. Run at multiple ε_org values. This introduces the multiphase transport enhancement (and ionic penalty). Compare polarisation curves and FE_ADPN to Stage 1 to quantify the D_mix effect.

Sub-steps: for each ε_org, warm-start from Stage 1 solution at V = −1.0 V → Newton continuation sweep over V → cache → plots → **STOP.**

**Stage 3: Full 3D sweep.**

Runs Stage 2 across the complete (ε_org, δ) grid. 30 Newton continuation sweeps. Generate performance map and all comparison plots → **STOP.**

---

## 13. Physicality Checks

| Check | Expected |
|-------|----------|
| φ_l(x) profile | A few mV variation; larger at high j |
| Σ z_i J_i at each face | = 0 to machine precision |
| Electroneutrality | \|Σ z_i c_i\| < 10⁻⁸ mol m⁻³ |
| c_AN(x) at high j | Depletes toward electrode |
| Buffer residuals at OCV | < 10⁻⁶ mol m⁻³ s⁻¹ |
| FE_ADPN vs ε_org | Should increase (enhanced AN transport) |
| D_AN,mix vs ε_org | Increases (organic pathway) |
| D_OH,mix vs ε_org | Decreases (volume exclusion) |
| No R_PT residual | Confirm R_PT is absent from residual |
| Bulk pH at x = δ | ~13.0 (matches solve_phosphate_equilibrium) |

---

## 14. Module Structure

```
an_ehd/
├── params.jl           # Constants (T, F, R, K_w, K_a2, K_a3, D_aq, D_org), Ref-backed kinetics
├── mesh.jl             # make_mesh(N, delta; stretch) — stretch-factor graded mesh
├── diffusivity.jl      # D_mix(i, eps_org), with m_i upgrade hook
├── chemistry.jl        # solve_phosphate_equilibrium, buffer_sources!, make_initial_guess
├── kinetics.jl         # j_ADPN, j_PN, j_HER
├── transport.jl        # sg_flux
├── assembly.jl         # full_residual! (no R_PT); integrated FV form
├── solver.jl           # newton_solve! (banded FD Jacobian, sparse LU); newton_continuation
├── run_stage1.jl       # Single-phase reference → STOP
├── run_stage2.jl       # D_mix(ε_org) sweep → STOP
├── run_stage3.jl       # Full 3D sweep → STOP
├── plot_results.py
└── output/cache/
```

Note: `solve_phosphate_equilibrium` and `make_initial_guess` live in `chemistry.jl` alongside `buffer_sources!`. `solver.jl` uses `SparseArrays.jl` + `KLU.jl` (or `LinearAlgebra.lu` with banded storage) — **no dense Jacobian**. No `phase_transfer.jl` module — R_PT is absent from this formulation.

---

## 15. Common Pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Using D_aq everywhere | ε_org has no effect | Use D_mix(ε_org) |
| D_org ≠ 0 for ions | Unphysical ionic transport | Set D_org = 0 for all ions |
| Adding ε_aq to R_buf | Double-counting | No prefactor (c_i is aqueous) |
| Including R_PT | Inconsistent with local equilibrium | Remove R_PT entirely |
| Wrong AN order | FE insensitive to j | c² for ADPN, c¹ for PN |
| Cold start (no equilibrium IC) | Newton diverges immediately | Use solve_phosphate_equilibrium + make_initial_guess |
| SG overflow | exp crash | Clamp α to [−700, 700] |
| Arithmetic D too weak | FE enhancement too small | Upgrade to m_i-corrected D_eff (§4.2) |
| Divided FV form on graded mesh | Jacobian ill-conditioned (entries O(1/dx_min)) | Use integrated form J_L − J_R + S·dx (§10.5) |
| Missing φ_l gauge fix | Jacobian rank-deficient, Newton stalls | Set φ_l(δ) = 0 (Dirichlet, §7.2) |
| j₀ in mA/cm² passed directly | Rates off by factor of 10 | Convert: A/m² = mA/cm² × 10 |
| Dense FD Jacobian | ~60 ms/step, full sweep >10 min | Banded FD with 19-column grouping + sparse LU (§10.2) |
| Fixed absolute dx_min/dx_max in mesh | Domain length wrong across δ sweep | Use stretch-factor parameterisation via make_mesh (§10.1) |
| φ_l step clamp at 0.1 V | Over-large potential steps destabilise c_i DOFs | Clamp φ_l steps to 0.015 V (§10.3) |
| Implementing PAC before testing Newton cont. | Unnecessary complexity; no fold bifurcations expected | Use simple Newton continuation first (§10.4) |

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

**Bloomquist et al.** (CEJ 2026, 528, 172125): j = 70–300 mA cm⁻², ε_org = 0.05–0.30. FE_ADPN = 73–76% at j > 200 mA cm⁻². Bubble convection dominates.

**Mathison et al.** (JACS 2025, 147, 4296): Mechanism — radical coupling for ADPN, proton transfer for PN (strong KIE).

**Suwanvaipattana et al.** (J. Cleaner Prod. 2017, 142, 1296): D, m, d_p values.

**Huang et al.** (CEJ 2020, 382, 123006): Langmuir-adsorption kinetics on Pb.

**Costentin & Savéant** (J. Electroanal. Chem. 564, 2004, 99): ΔG⁰ = −1.84 eV for protonated radical coupling.

---

## 18. Fitting Strategy

1. **ε_org = 0, δ = 50 μm.** Fit kinetic parameters (j₀, α_c) to single-phase FE data.
2. **Activate D_mix at ε_org = 0.15.** Check whether FE_ADPN increases. If not: the arithmetic D_mix is too weak → upgrade to m_i-corrected D_eff.
3. **δ sweep at ε_org = 0.15.** Match j_total trends to Bloomquist data.
4. **ε_org sweep at fixed δ.** Match FE vs organic loading curve.
5. **Full 3D sweep.** Generate performance map.

| Metric | Target |
|--------|--------|
| FE_ADPN peak (ε_org = 0) | 50–60% |
| FE_ADPN peak (ε_org = 0.15) | 73–80% |
| FE_ADPN at j = 200 mA cm⁻² | > 70% |
| FE enhancement from ε_org = 0 → 0.15 | +20–30 pp |

If the arithmetic D_mix cannot produce this enhancement, the m_i correction (§4.2) is needed.

---

## 19. Required Plots

| Panel | Content | Axes |
|-------|---------|------|
| (a) | j_ADPN, j_PN, j_HER vs V | −V vs SHE / mA cm⁻² |
| (b) | FE_ADPN, FE_PN, FE_HER vs V | −V vs SHE / % |
| (c) | FE_ADPN vs j at multiple ε_org | j / FE — key Bloomquist comparison |
| (d) | FE_ADPN vs ε_org at fixed j | ε_org / FE |
| (e) | c_AN(0)/c_AN,bulk vs j | j / depletion ratio |
| (f) | φ_l(0) vs j at multiple ε_org | j / mV (ohmic penalty) |
| (g) | D_AN,mix and D_OH,mix vs ε_org | ε_org / D |
| (h) | Production rate vs ε_org | ε_org / kg cm⁻² h⁻¹ |

---

*References: Bloomquist et al. CEJ 2026; Corpus et al. Joule 2023; Weng, Bell & Weber PCCP 2018; Huang et al. CEJ 2020; Suwanvaipattana et al. J. Cleaner Prod. 2017; Mathison et al. JACS 2025; Costentin & Savéant J. Electroanal. Chem. 2004; Lasia J. Electroanal. Chem. 1995.*
