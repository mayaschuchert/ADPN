# ADPN EHD Implementation Guide — Changelog v4 → v5

Revision date: 2026-04-21

This changelog lists every substantive difference between `ADPN_EHD_Implementation_Guide_v4.md` and `ADPN_EHD_Implementation_Guide_v5.md`. Changes are grouped by topic. Each entry notes the section(s) affected and flags whether the change is a **physics** change (model behaviour differs), a **numerical** change (solver/infrastructure), a **documentation** clarification, or a **sweep-range** change.

---

## 1. Bulk acrylonitrile concentration — Convention A per total volume (physics)

**Affected sections:** §2 (intro), §4 (table), §7.2, §8 (sweep range), §9.3 (electrolyte composition), §12 (Stage 1), §13 (physicality), §15 (pitfalls), §18 (fitting), §19 (plots).

### 1.1 Formula corrected

- **v4:** `c_AN,bulk = C_AN,total / (1 − ε_org)` for ε_org < 0.09, with `C_AN,total` implicitly defined. The v4 code implemented this as `C_AN_SAT / (1 − ε_org)` — which unconditionally pins the numerator at the saturation concentration. That is physically wrong: at ε_org = 0.05 it gives c_AN = 1379 mol/m³, which is *above* the aqueous saturation (a one-phase solution cannot exceed saturation).
- **v5:** `c_AN,bulk = ε_org · ρ_AN / M_AN` (Convention A — moles of AN per total solution volume) for ε_org < ε_sat, pinned at `C_AN_SAT = 1310 mol/m³` for ε_org ≥ ε_sat.

### 1.2 ε_sat derived physically

- **v4:** ε_sat = 0.09 (quoted, not derived).
- **v5:** ε_sat = C_AN_SAT / (ρ_AN/M_AN) = 1310 / 15191 ≈ **0.0862** (derived from the saturation continuity condition ε_sat · ρ_AN/M_AN = C_AN_SAT). v4's value 0.09 over-rounds this.

### 1.3 Constants added

Guide explicitly introduces AN physical constants in §7.2 and §9:

| Symbol | Value | Unit | Source |
|---|---|---|---|
| M_AN | 0.05306 | kg/mol | molar mass of C₃H₃N |
| ρ_AN | 806 | kg/m³ | neat liquid AN at 25 °C (Suwanvaipattana 2017) |
| m_AN | 11.59 | — | partition coefficient (Suwanvaipattana 2017) |
| ρ_AN/M_AN | **15 191** | mol/m³ | molar density of neat AN |
| C_AN_SAT | **1 310** | mol/m³ | = ρ_AN/(M_AN·m_AN), aqueous saturation |
| EPS_ORG_SAT | **0.0862** | — | = C_AN_SAT / (ρ_AN/M_AN) |

### 1.4 ε_org = 0 now truly non-physical

- **v4:** ε_org = 0 in Stage 1 was "single-phase reference, no organic droplets, D = D_aq." The AN bulk was still 1310 mol/m³ due to the wrong formula.
- **v5:** ε_org = 0 gives c_AN,bulk = 0 — no AN in the system, no reduction reactions possible. `log(c_AN) = −Inf` crashes the initial guess. ε_org = 0 is flagged as pathological throughout.

### 1.5 Sweep range shifted

- **v4:** `ε_org ∈ {0.00, 0.05, 0.09, 0.15, 0.25, 0.30}` (guide §8).
- **v5:** `ε_org ∈ {0.02, 0.05, 0.08, 0.15, 0.25, 0.30}` — spans both regimes with 0.02 as the minimum physical loading. Three points below ε_sat (single-phase), three above (two-phase). ε_org = 0 removed from the sweep list.

### 1.6 Stage 1 anchor moved

- **v4:** Stage 1 runs at ε_org = 0 (degenerate with the corrected physics).
- **v5:** Stage 1 runs at **ε_org = 0.02** — minimum physical single-phase reference. `c_AN,bulk ≈ 304 mol/m³`, D_mix = D_aq for all species.

---

## 2. Regime-aware effective diffusivity D_mix (physics)

**Affected sections:** §4, §13 (physicality), §15 (pitfalls), §19 (plot g).

### 2.1 Single-phase / two-phase split

- **v4:** Arithmetic mean applied universally: `D_i,mix = ε_org · D_i,org + (1 − ε_org) · D_i,aq` for all ε_org. Below solubility, this treated the "dissolved AN" as a parallel pathway that doesn't physically exist.
- **v5:** Regime-aware formula:

```
D_i,mix = D_i,aq                                          if ε_org < ε_sat
        = ε_org · D_i,org + (1 − ε_org) · D_i,aq          if ε_org ≥ ε_sat
```

Below ε_sat there are no organic droplets — the system is a single aqueous phase with dissolved AN, and transport proceeds through water at D_i,aq. Only above ε_sat does the volume-averaged mixing rule apply.

### 2.2 Regime transition now a physicality check

**v5:** New physicality check in §13: "FE_ADPN vs ε_org — peaks around ε_sat ≈ 0.086 ↑ in two-phase regime." Plot (g) in §19 must show the step transition in D_AN,mix and D_OH,mix at ε_sat.

---

## 3. Phosphate buffer chemistry — OH⁻-pathway (physics)

**Affected sections:** §6 (entirely rewritten), §13 (physicality check for pH profile added), §15 (new pitfall row).

### 3.1 Reactions reformulated

- **v4:** H⁺-pathway:
  - R1: H₂O ⇌ H⁺ + OH⁻
  - R2: H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻
  - R3: HPO₄²⁻ ⇌ H⁺ + PO₄³⁻

- **v5:** OH⁻-pathway (kinetically more realistic at pH > 10):
  - R1: H₂O ⇌ H⁺ + OH⁻ (unchanged stoichiometry, new constants)
  - R2: **OH⁻ + H₂PO₄⁻ ⇌ HPO₄²⁻ + H₂O** (new)
  - R3: **OH⁻ + HPO₄²⁻ ⇌ PO₄³⁻ + H₂O** (new)

### 3.2 Source stoichiometry change

| Species | v4 (H⁺-pathway) | **v5 (OH⁻-pathway)** |
|---|---|---|
| H⁺ | r₁ + r₂ + r₃ | **r₁** |
| OH⁻ | r₁ | **r₁ − r₂ − r₃** |
| H₂PO₄⁻ | −r₂ | −r₂ (unchanged) |
| HPO₄²⁻ | r₂ − r₃ | r₂ − r₃ (unchanged) |
| PO₄³⁻ | r₃ | r₃ (unchanged) |

Charge conservation Σ zᵢ·Sᵢ = 0 holds in both formulations.

### 3.3 Rate constants updated

| Constant | v4 | **v5** | Unit |
|---|---|---|---|
| Ka2 | 6.2 × 10⁻⁵ | **6.3 × 10⁻⁵** | mol/m³ |
| Ka3 | 4.8 × 10⁻¹⁰ | **4.5 × 10⁻¹⁰** | mol/m³ |
| k₁,f | 2.4 × 10⁻⁵ (labeled "s⁻¹" — dimensionally ambiguous) | **1.4 mol/(m³·s)** (Eigen–De Maeyer 1955) | mol/(m³·s) |
| k₁,r | = k₁,f / K_w ≈ 2 400 | = k₁,f / K_w = **1.4 × 10⁸** | m³/(mol·s) |
| k₂,f | 10⁶ s⁻¹ (first-order forward) | **10⁵ m³/(mol·s)** (second-order forward) | m³/(mol·s) |
| k₂,r | derived | derived ≈ **15.87** | s⁻¹ |
| k₃,f | 10² s⁻¹ | **2 × 10³ m³/(mol·s)** | m³/(mol·s) |
| k₃,r | derived | derived ≈ **4.44 × 10⁴** | s⁻¹ |

### 3.4 Bulk pH shifts slightly

Consequence of the Ka2, Ka3 updates:

- **v4:** bulk pH ≈ **13.016** at 0.5 M Na₃PO₄ + 0.02 M TBA-OH.
- **v5:** bulk pH ≈ **13.027** (shift of +0.011).

Negligible for downstream predictions but worth noting for continuity of results across versions.

### 3.5 Surface pH anomaly resolved

**v4 with H⁺-pathway:** at high current (V = −2.5 V, j ~ 10 A/cm²), the model produced *pH at surface < pH at bulk* — counter-intuitive because OH⁻ is produced at the electrode and should drive pH up. Root cause: k₁,f = 2.4 × 10⁻⁵ s⁻¹ is too slow for water autoprotolysis to re-equilibrate K_w locally when c_OH jumps 2000× at the surface.

**v5 with OH⁻-pathway and Eigen–De Maeyer k₁,f:** at the same conditions, surface pH ≈ 16 and bulk pH ≈ 13 — the physically-expected direction. The fast OH⁻-pathway phosphate kinetics and the 60× larger water rate jointly keep all equilibria within ~machine precision of their K_eq values at moderate currents.

New pitfall row in §15:
> H⁺-pathway buffer at high pH → Unphysical k_r ~ 10¹⁰ m³/(mol·s) → **Fix:** OH⁻-pathway stoichiometry (§6)

---

## 4. Numerical method updates

### 4.1 Newton solver: direct `(J+λI) du = −F` (numerical)

**Affected sections:** §10.4, §15 (two new pitfalls).

- **v4:** Damped Newton `J du = −F` with row scaling and step clamping; constant `λ = 0` (pure Newton).
- **v5:** Direct `(J + λI) du = −F` with fixed Tikhonov `λ = 10⁻¹⁰` (keeps near-singular Jacobians solvable without distorting the Newton direction). **Strict L2 descent** with 10-step backtracking line search. Same per-DOF step clamps (log-conc ≤ 5.0, φ_ℓ ≤ 0.015 V).

Documented in §15:
> Classical LM `(JᵀJ+λI)du=−JᵀF` with λ₀=1 → Picks wrong solution branch → **Fix:** Use direct `(J+λI)du=−F` with λ=10⁻¹⁰ (§10.4)
>
> Merit slack 1.10 at tight tol → Accepts numerical-noise as progress → **Fix:** Strict L2 descent at tol ≤ 10⁻⁵

The v5 code retains the classical LM path only as a diagnostic — the production solver uses direct Newton.

### 4.2 Convergence tolerance loosened (numerical)

- **v4:** `‖F‖_∞ < 10⁻⁴`.
- **v5:** `‖F‖_∞ < 10⁻⁵` as the default. Charge conservation `|Σ z_i N_i|_interior < 10⁻¹⁰ mol/m²/s` is the real correctness gate (physicality §13) — tightening `tol` below 10⁻⁵ wastes iterations chasing unreachable precision in near-singular regions (cond(J) ~ 10¹⁶).

### 4.3 Jacobian half-bandwidth corrected: 9 → 17 (numerical, critical)

**Affected sections:** §10.2, §15.

- **v4:** stated `b = 9` (the block size, not the bandwidth), with FD column-grouping stride 2b + 1 = 19.
- **v5:** **b = 2·B − 1 = 17** for block-tridiagonal coupling with block size B = 9 in cell-major layout. Column-grouping stride 2b + 1 = **35**.

⚠️ **With v4's `b = 9` the banded FD Jacobian was silently incorrect** — unrelated columns at stride 19 aliased onto the same rows, producing an exactly-singular Jacobian. This was diagnosed by comparing the banded FD to a dense FD Jacobian at runtime; they disagreed only with `b = 9`, and agreed to machine precision with `b = 17`.

New pitfall row:
> Jacobian halfbandwidth b = 9 → Banded FD aliases columns, exactly singular → **Fix:** Use b = 2·block_size − 1 = 17 (§10.2)

### 4.4 ForwardDiff AD Jacobian option (numerical)

**New sections:** §10.3, §10.8 (type-generic residual chain).

- **v5:** Added `jacobian_mode = :ad` option to `newton_solve!` that uses `ForwardDiff.jacobian!` with `JacobianConfig` (chunk size 12). Requires the entire residual chain to be type-generic (`AbstractVector{<:Real}`, `::Real` for scalars, `zeros(eltype(u), ...)` for work arrays).

Both modes (`:fd` and `:ad`) are kept. FD is ~4× faster per Jacobian and adequate in most regions. AD is exact to machine precision and more robust near rank-deficiency (CO2R-guide-style Brady-et-al. effect). In direct-Newton mode with the Taylor-smoothed Bernoulli (§4.5), the two modes give identical converged solutions.

### 4.5 Taylor-smoothed Bernoulli branch in SG flux (numerical, critical for AD)

**Affected sections:** §3.3, §15 (new pitfall).

- **v4:** Hard branch in `sg_flux`:
  ```julia
  if abs(alpha) < 1e-10
      return -D * (c_R - c_L) / dx
  end
  ```
- **v5:** Smooth Taylor expansion for `|α| < 0.01`:
  ```julia
  if abs(alpha) < 0.01
      B_pos = 1 − α/2 + α²/12
      B_neg = 1 + α/2 + α²/12
  else
      # full Bernoulli
  end
  ```

Required because ForwardDiff AD, evaluated at the α = 0 equilibrium state, goes down the short-circuit branch and computes `dJ/dφ = 0` (wrong) — the hard branch has no φ-dependence in that region. FD perturbation trips into the other branch and gets the correct derivative. With the Taylor form, both FD and AD compute identical (correct) derivatives through α = 0.

New pitfall row:
> SG hard branch at α=0 → AD derivative spuriously 0 → **Fix:** Taylor expansion for |α| < 0.01

### 4.6 Type-generic residual chain (numerical)

**Affected sections:** §10.8.

- **v4:** `full_residual!(res, u, ...)` with `res::AbstractVector`, `u::AbstractVector`. Inner function signatures and work-array allocations (`zeros(8, N)`) were `Float64`-typed.
- **v5:** `full_residual!(res::AbstractVector{T}, u::AbstractVector{T}, ...) where {T<:Real}`. Work arrays use `zeros(T, 8, N)`. `sg_flux` and `tafel_currents` signatures relaxed from `::Float64` to `::Real` for scalar arguments. Enables ForwardDiff.Dual numbers to propagate.

### 4.7 Optional log-j continuation (numerical)

**Affected sections:** §10.5.

- **v4:** Simple Newton continuation in V only, with adaptive ds.
- **v5:** Same V-continuation as the default (reaches V = −2.5 V reliably in single-phase Stage 1). Adds an optional `newton_continuation_logj` that parameterises the path by log₁₀(j_total) via augmented Newton with V as an extra DOF. Useful for dense sampling in the exponential-Tafel regime. Not enabled by default since the V floor is set by model validity, not step-size control.

### 4.8 Minor adaptive-continuation refinements

**Affected sections:** §10.5.

- **v4:** On success in ≤ 5 iter → ds × 1.5; on failure → ds / 2.
- **v5:**
  - Graduated growth: ≤ 4 iter → ×1.4, ≤ 10 iter → ×1.1, ≤ 20 iter → ×1.0, > 20 iter → ×0.7.
  - Aggressive shrink on failure: ×0.3.
  - `max_total_fail = 200` safety cap.
  - `ds_min` tightened from 0.005 V to **10⁻⁴ V**.

### 4.9 Roots.jl dependency removed (numerical, Windows-specific)

**Affected sections:** §6.4, §14.

- **v4:** `solve_phosphate_equilibrium` uses `Roots.find_zero(..., Bisection())`.
- **v5:** Inline bisection routine (10 lines). Reason: Windows Defender Application Control blocks `Roots.jl`'s precompiled DLL cache on some systems (error `0xc0e90002`). The inline bisection is equivalent for the simple 1D charge-balance problem and has no external dependency.

---

## 5. Plot layout updates (documentation)

**Affected sections:** §19.

### 5.1 Profile figures: 2×2 → 2×3 with pH panel

- **v4:** 2×2 grid: H⁺/OH⁻, phosphates, AN/ADPN/PN, φ_ℓ.
- **v5:** 2×3 grid:

| (0,0) H⁺ and OH⁻ (log y) | (0,1) Phosphate speciation | (0,2) **pH** (new) |
| (1,0) AN / ADPN / PN | (1,1) φ_ℓ | (1,2) reserved |

The pH panel shows `pH(x) = −log₁₀(c_H/1000)` with a dashed reference line at the bulk pH. This is the most direct way to see whether the electrode-adjacent region is more or less basic than bulk — a sensitive diagnostic of buffer kinetics (the surface-pH anomaly that motivated the OH⁻-pathway change in §3).

### 5.2 Polarisation panel split: (a) → (a1) log + (a2) linear

- **v4:** Panel (a): j_r vs V on log axis only.
- **v5:** Panel **(a1)** log-scale polarisation + **Panel (a2)** linear-scale polarisation. The linear view exposes the exponential Tafel "takeoff" structure that's compressed to a straight line on the log axis.

### 5.3 D_AN,mix / D_OH,mix regime-transition plot

- **v4:** Panel (g) showed D vs ε_org as a smooth curve (arithmetic mean).
- **v5:** Panel (g) shows the **step transition at ε_sat ≈ 0.0862** — flat at D_aq below, arithmetic mean above. This is now a required diagnostic of the regime-aware D_mix.

---

## 6. Sweep and stage structure (documentation + physics)

### 6.1 Three-parameter sweep table

- **v4:** ε_org = {0.00, 0.05, 0.09, 0.15, 0.25, 0.30}.
- **v5:** ε_org = **{0.02, 0.05, 0.08, 0.15, 0.25, 0.30}**. Single-phase points below ε_sat and two-phase points above, spanning both regimes. ε_org = 0 explicitly banned.

### 6.2 Stage 1 anchor

- **v4:** `STAGE1_EPS = 0.00` — validates numerics, no AN effect.
- **v5:** `STAGE1_EPS = 0.02` — minimum-physical single-phase reference. Still D_i,mix = D_i,aq everywhere (single-phase regime).

### 6.3 Fitting strategy anchor

- **v4:** Step 1 of §18 fit: "ε_org = 0, δ = 50 μm."
- **v5:** Step 1: "**ε_org = 0.02, δ = 50 μm.**"

### 6.4 Experimental targets (Bloomquist comparison)

- **v4:** "FE enhancement from ε_org = 0 → 0.15 = +20–30 pp."
- **v5:** "FE enhancement from **ε_org = 0.05 → 0.15** = +20–30 pp" — both points now physical, both in distinct regimes (single-phase vs two-phase).

---

## 7. New or expanded physicality checks (§13)

Added / modified:

- **"pH(x): Rises monotonically from bulk to surface at high j (OH⁻-pathway)"** — directly tests whether the new buffer chemistry behaves correctly. Replaces the implicit v4 expectation.
- **"FE_ADPN vs ε_org: Peaks around ε_sat ≈ 0.086 ↑ in two-phase regime"** — explicit regime-transition diagnostic.
- **"D_AN,mix vs ε_org: Flat (single-phase) then rises (two-phase)"** and the analogous OH⁻ check — explicitly reference the regime split.
- Bulk pH expected value: "≈ 13.016" → "**≈ 13.03**" (Ka update).

---

## 8. Module structure (documentation)

**Affected sections:** §14.

- **v5** adds comments indicating:
  - `params.jl` exports `MOLAR_DENSITY_AN`, `C_AN_SAT`, `EPS_ORG_SAT`
  - `chemistry.jl` uses the OH⁻-pathway buffer and Convention A `c_AN_bulk`
  - `transport.jl` uses Taylor-smoothed Bernoulli for |α| < 0.01
  - `assembly.jl` is type-generic `AbstractVector{T<:Real}`
  - `solver.jl` exports both `:fd` and `:ad` Jacobian modes
  - `solve_phosphate_equilibrium` uses inline bisection (no `Roots.jl`)

---

## 9. New pitfall rows in §15

(All new in v5; highlights the main places a naïve implementation will fail.)

| Pitfall | Fix |
|---|---|
| ε_org = 0 in sweep | Sweep from ε_org ≥ 0.02 |
| Single-phase with mixed D | D_mix = D_aq for ε_org < ε_sat |
| H⁺-pathway buffer at high pH | OH⁻-pathway stoichiometry (§6) |
| SG hard branch at α=0 | Taylor expansion for \|α\| < 0.01 |
| Jacobian halfbandwidth b = 9 | Use b = 2·block_size − 1 = 17 |
| Classical LM `(JᵀJ+λI)du=−JᵀF` with λ₀=1 | Use direct `(J+λI)du=−F` with λ=10⁻¹⁰ |
| Merit slack 1.10 at tight tol | Strict L2 descent at tol ≤ 10⁻⁵ |
| Non-type-generic residual | `AbstractVector{T<:Real}`; `zeros(eltype(u),...)` |

---

## 10. New references cited

**§17.**

- **Eigen & De Maeyer** (Z. Elektrochem. 1955): water autoprotolysis rate k₁,f = 1.4 × 10⁻³ M/s (= 1.4 mol/(m³·s)). Anchors the v5 k₁,f constant.

---

## Net impact summary

| Category | v4 → v5 change |
|---|---|
| Physics — AN bulk | Corrected (Convention A per total volume); ε_org = 0 forbidden |
| Physics — D_mix | Regime-aware (single-phase = D_aq, two-phase = arithmetic mean) |
| Physics — buffer | OH⁻-pathway with Eigen-De Maeyer water rates |
| Physics — sweep | ε_org ∈ [0.02, 0.30], never 0 |
| Numerics — Jacobian | Bandwidth corrected 9 → 17; AD option via ForwardDiff |
| Numerics — SG flux | Taylor-smooth through α = 0 |
| Numerics — Newton | Direct `(J+λI)du=−F` + strict L2 descent + tol = 10⁻⁵ |
| Numerics — continuation | Graduated adaptive ds; optional log-j variant |
| Numerics — deps | `Roots.jl` removed (inline bisection) |
| Documentation — plots | 2×3 profile grid with pH panel; a1/a2 polarisation split |
| Documentation — stages | Stage 1 moved from ε_org = 0 to ε_org = 0.02 |

Stage 1 at ε_org = 0.02 with the v5 pipeline reaches V = −2.5 V reliably in 15 continuation points, with all physicality checks passing to the expected precision (charge conservation at 10⁻¹¹ mol/m²/s, buffer @ bulk at 10⁻⁹, pH monotone increase toward surface under HER).
