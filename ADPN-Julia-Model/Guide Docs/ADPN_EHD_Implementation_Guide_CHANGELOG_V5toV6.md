# ADPN EHD Implementation Guide — CHANGELOG v5 → v6

**Status:** v6 written 2026-04-27. v5 remains the canonical reference for the *internal* model (NP transport, Tafel kinetics, OH⁻-pathway buffer, regime-aware D_mix, Newton solver). v6 adds the *external* coupling needed to compare the model with Bloomquist et al. (CEJ 2026) experimental data. Nothing in v5's solver, residual, or chemistry has been altered.

---

## 1. Why v6 exists

After Stages 2 and 2m the model peaks at FE_ADPN ≈ 38–49% (m_i-corrected D_eff) versus Bloomquist's 73–76% headline. The dominant lever is kinetics, not transport — but kinetics cannot be fit until the model and the experiment share an axis. Bloomquist controls *current density at fixed cell voltage*, while the v5 model solves *at fixed cathode potential vs SHE*. v6's job is to bridge these two pictures with the smallest physically defensible scaffolding so that the six Tafel parameters can finally be fit against 162 rows of real data.

---

## 2. Sections added or rewritten

| § | Title | Status |
|---|---|---|
| 17 | Cell-Voltage Decomposition | NEW |
| 18 | Hydrodynamics: Flow → δ Mapping | NEW |
| 19 | Experimental Data | rewritten — Bloomquist CSV schema added |
| 20 | Fitting Strategy | rewritten — kinetics-only, transport frozen |
| 21 | Required Plots | extended — fit-validation panels (i)–(o) |
| 12 | Implementation Stages | Stage 3 reframed as "warm-start cache + default-kinetics forward sweep on Core"; Stage 4 split into Stage 4a (Core fit) + Stage 4b (Extended/Holdout forward apply, no re-fit); Stage 4c (joint refinement) explicitly deferred to v7 |
| 14 | Module Structure | `cell_voltage.jl`, `hydrodynamics.jl`, `fixed_j_solver.jl`, `fit_kinetics.jl`, `data/`, `plot_fit.py`, `run_stage4.jl` added |
| 15 | Pitfalls | nine v6-specific pitfalls added |
| Header / Footer | scope sentence; v6 provenance line | updated |

Sections 1–16 (governing equations, kinetics, buffer chemistry, BCs, parameter tables, numerical methods, caching, physicality checks, potential referencing) are **unchanged** from v5.

---

## 3. Design decisions (and their explicit simplifications)

This is the part to read carefully — every choice below is a place where v6 punts something that the underlying physics actually wants.

### 3.1 Cell-voltage decomposition uses a single lumped V_CE

```
V_cell_meas = V_CE + |V_cathode_SHE| + j · [(gap − δ)/κ_eff + R_contact]
```

| What V_CE absorbs | Why we lumped it |
|---|---|
| E°_OER thermodynamic potential | Constant; no fit-leverage |
| Average η_anode at the j-range of interest | Tafel-form would add 2 params; data doesn't constrain them once R_contact is also fit |
| Anode-side mass-transport overpotential | Small for SS-OER in this regime; not separable from η_anode in two-electrode data |
| Anode-side double-layer / film effects | Not separately observable |
| **What it doesn't absorb:** the j-dependence of η_anode beyond a constant | Would need Tafel(η_a) → goes into v7 if residuals demand it |

Initial guess V_CE = 1.7 V (= 1.23 V thermodynamic OER + 0.45 V typical SS overpotential at ∼100 mA cm⁻²). Bounds [1.4, 2.1] V.

### 3.2 Electrolyte conductivity computed, not fit

```
κ_dilute = (F²/RT) · Σ z_i² · D_i,aq · c_i,bulk    (§17.2)
κ_eff    = κ_dilute · (1 − ε_org)^1.5              (Bruggeman)
```

Fed by the model's own bulk equilibrium (§6.4): H⁺, OH⁻, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻, Na⁺ at pH 13.03. Predicts κ ≈ 6.6 S/m at ε_org = 0 dropping to ≈ 3.8 S/m at ε_org = 0.30 — same order of magnitude as the empirical 5–10 S/m range in CONTEXT_TRANSFER §7.

| What we punt | Reason |
|---|---|
| Concentrated-solution corrections (activity coefficients, ion association) | Phosphate at 0.5 M is borderline; dilute theory off by maybe 20% — small relative to fit residuals |
| Migration enhancement of κ from gradients | Already inside `V_cathode_SHE` via φ_l(0); the (gap − δ) bulk term assumes well-mixed |
| Bubble void Bruggeman `(1 − ε_gas)^1.5` | **Deferred to v7** — the dominant missing physics for the 0.25 mm gap |
| TBA⁺ vs Na⁺ mobility distinction | TBA is lumped into Na⁺ (v5 §2 convention); modest κ overestimate, ~1 % |

### 3.3 Lévêque correlation, no bubble enhancement

```
Sh = 1.85 · (Re·Sc·d_h/L)^(1/3),   δ_lam = d_h / Sh   (§18)
```

W = 4 mm (channel), L = 0.16 m (= A/W), ν = 10⁻⁶ m²/s. δ ranges 60–310 μm across all 9 (gap × Q_total) blocks.

| What we punt | Reason |
|---|---|
| K_δ geometric correction (CONTEXT_TRANSFER proposed [0.3, 3.0]) | Set K_δ = 1.0 in v6; promote to fit param in v7 if residuals show systematic Q-bias |
| Bubble-induced convection enhancement of δ | **Deferred to v7** — the headline Bloomquist finding ("bubbles dominate"); v6 will systematically *underestimate* mass transport at high j |
| Two-phase BL (organic droplet contribution) | Lumped into the regime-aware D_mix (§4), not into δ |
| Serpentine entrance/exit effects | Subsumed into single L; reasonable when L/d_h ≫ 1 |
| Channel width W = 4 mm fixed | Would only matter if we got systematic gap-dependent FE bias *after* kinetic fit; tier-2 fit param if needed |

The Lévêque expression is valid for laminar flow with `Re·Sc·d_h/L ≫ 1`. All Bloomquist rows have Re < 10 and Sc ≈ 435, so the developing-BL assumption holds. We are *well* below the laminar-turbulent transition.

### 3.4 Bubble physics deferred wholesale

This is the largest single simplification. v6 contains *no* bubble term anywhere — not in δ, not in κ_eff, not in V_CE. Justification: the Bloomquist headline result is that bubble-induced convection ≈ raises FE_ADN by 10–25 pp at high j in droplet-like regimes. v6 will reproduce the *shape* of FE vs (j, ε_org) but undershoot magnitude on the 0.25 mm gap holdout. This is intentional: better to have a kinetics fit that is interpretable on transparent transport assumptions than one that buries kinetic error inside a tunable bubble parameter.

**v7 entry points (already stubbed in v6 docstrings):**
- `delta_actual = delta_leveque · f_bubble(j, gap, Q)` in `hydrodynamics.jl`
- `κ_eff *= (1 − ε_gas)^1.5` in `cell_voltage.jl`
- One scalar fit param `c_bubble` (Vogt-style) or two (`α_b`, `β_b` in `f_bubble = 1 + α_b · (j/j_ref)^β_b`)

### 3.5 TCH (tricyanohexane) excluded from the fit

v5 has 8 species; Bloomquist measures 4 products (ADN, TCH, PN, H₂). v6 fits residuals on (FE_ADN, FE_PN) only — `FE_HER ≈ 100% − FE_ADN − FE_TCH − FE_PN` is reconstructed for diagnostic but not fit. TCH FE is typically 5–17%, so the FE_HER residual after summing has a mean of ~10% systematic offset. v7 candidate: add TCH as species 9 with `j_TCH ∝ c_AN³` Tafel.

### 3.6 Three-tier row selection: Core / Extended / Holdout

v6 does *not* fit on all 162 rows. Including rows where the model is structurally wrong (missing bubble physics, missing TCH species, AN-starved degeneracy) would pull kinetic parameters toward unphysical values. Three concentric subsets, ordered by trust in v6 physics:

| Subset | Filter | Rows | Used for |
|---|---|---|---|
| **Core** | gap ∈ {0.5, 1.0} mm AND j ≤ 190 mA cm⁻² AND ε_org ≥ 0.04 | ≈60 | Stage 4a fit |
| **Extended** | gap ∈ {0.5, 1.0} mm AND ε_org ≥ 0.04 (any j) | ≈96 | Stage 4b forward apply, no re-fit |
| **Full holdout** | gap = 0.25 mm AND ε_org ≥ 0.04 | ≈48 | Stage 4b forward apply, untouched during fitting |

| Filter | Rows dropped | Justification |
|---|---|---|
| ε_org < 0.04 (all subsets) | 18 of 162 (11%) | AN-starved regime; FE_ADN ≈ 0 experimentally; physically degenerate, adds noise without information |
| j > 190 mA cm⁻² (Core only) | excluded from Core, included in Extended | High j → bubble convection dominates; v6's Lévêque δ has no bubble enhancement → systematic FE_ADN under-prediction |
| gap = 0.25 mm (Core + Extended) | 54 of 162 (33%) | Bubbles span the gap; void blocking dominates ohmic drop *and* mass transport; no v6 model term captures this |

Stage 4a fits **only on Core**. Stage 4b forward-applies the converged params to Extended and Full holdout *without re-fitting*. The two forward-apply RMSE deltas measure how badly bubble physics is missing:

- Core RMSE → "is v6 self-consistent on its own valid envelope?"
- Extended RMSE − Core RMSE → "how much does high-j bubble convection cost us in mass transport?"
- Holdout RMSE − Extended RMSE → "how much does small-gap bubble void blocking cost us in ohmic drop?"

Each gate has a numeric threshold in §20.4. Crossing them triggers specific v7 work, not v6 fit re-tuning.

> **Why fitting on all 162 rows would be worse, even though more data is "more information."** With a 6-param fit and 324 residuals, statistical power is not the constraint. Structural model error *is*. An optimiser presented with 54 0.25 mm gap rows that require bubble physics will move (j₀,r, α_c,r) to absorb the missing bubble correction — yielding kinetic parameters that fit 0.25 mm well at the cost of being wrong everywhere else. Sequestering the structurally-wrong rows into holdout keeps the Core fit interpretable.

### 3.7 Fixed-j solver replaces fixed-V continuation for fitting

Bloomquist data is constant-current. v6 introduces `fixed_j_solver.jl` which bisects on V vs SHE around the v5 Newton-continuation root such that `Σ j_r(V) = j_target`. Bracket [−2.5, −0.8] V. The v5 Newton solver is reused unchanged inside the bisection inner loop.

> Implementation note: warm-start V from the previous (gap, Q, ε_org, j) row's converged V to keep Newton's basin of attraction; for the first row of each (gap, Q, ε_org) block, run a quick V continuation from −1.0 V to bracket the root.

### 3.8 V_CE and R_contact frozen — only six kinetic params fit

v6 does **not** fit V_CE and R_contact. Reason: the Bloomquist SI tables omit per-row V_cell. V_cell can be back-derived from the energy productivity column `EP_ADN` but the reverse map compounds measurement noise from PR_ADN, j, and EP_ADN — fitting against this back-derived quantity would *reduce* the trustworthiness of V_CE and R_contact.

Two cleaner v7 paths once bubble physics lands:

1. Re-acquire per-row V_cell directly from a future Bloomquist data refresh (raw data exists).
2. Couple V_CE and R_contact to ε_gas(j) and η_anode(j), then jointly fit (kinetics + voltage + bubble) once the bubble model is in.

| Parameter | Fit in v6? | Source if not fit |
|---|---|---|
| j₀,1, j₀,2, j₀,3 | ✅ | — |
| α_c,1, α_c,2, α_c,3 | ✅ | — |
| V_CE | ❌ (frozen) | 1.7 V vs SHE — defer to v7 |
| R_contact | ❌ (frozen) | 1×10⁻⁴ Ω·m² — defer to v7 |
| E°_1, E°_2 (ADPN/PN onset) | ❌ | Mathison JACS 2025 → −1.30 V vs SHE |
| E°_3 (HER onset) | ❌ | Nernst at pH = 14 → −0.83 V vs SHE |
| D_i,aq, D_i,org | ❌ | Suwanvaipattana 2017, CRC |
| m_i (partition) | ❌ | Suwanvaipattana 2017 (used only if §4.2 m_i upgrade activated) |
| K_w, K_a2, K_a3, k_{1..3,f}, k_{1..3,r} | ❌ | Eigen 1955 + thermodynamic K_eq |
| C_AN_SAT, EPS_ORG_SAT | ❌ | Derived from ρ_AN, M_AN, m_AN |
| W (channel width) | ❌ | Fixed 4 mm |
| L (channel length) | ❌ | A / W = 16 cm |
| ν (kinematic visc.) | ❌ | 10⁻⁶ m²/s, water-like |
| σ_AN-water | ❌ | 10.5 mN/m (Girifalco–Good, SI) |
| ρ_aq, ρ_org | ❌ | 1000, 810 kg/m³ |

**Total fit dimension: 6.** Core subset has ≈60 rows × 2 residual species = ≈120 residuals. 20× overdetermined — comfortably enough for a well-posed LM fit.

---

## 4. New artefacts

### 4.1 Code modules (implemented in v6 — see §14 of guide for full tree)

```
an_ehd/
├── kinetics.jl          — patched: KIN_OVERRIDE Ref allows Stage 4 to vary
│                          (j₀, α_c) without mutating Params. Default = nothing,
│                          so Stages 1/2/2m/3 are byte-identical to v5.
├── hydrodynamics.jl     — NEW: d_hydraulic, v_super, delta_leveque,
│                          weber_numbers, reynolds, schmidt, sherwood_leveque,
│                          ml_min_to_m3_s
├── cell_voltage.jl      — NEW: kappa_dilute, kappa_eff, R_series,
│                          V_cell_predicted, V_cathode_target
├── fixed_j_solver.jl    — NEW: solve_at_j (bisects V vs SHE around the v5
│                          Newton solver, warm-started, optional KIN_OVERRIDE
│                          push/restore)
├── fit_kinetics.jl      — NEW: BloomquistRow, FitContext, build_context,
│                          select_core/extended/holdout, residuals!, loss,
│                          theta_to_physical, lm_fit (pure-Julia LM, no deps)
├── run_stage4.jl        — NEW: load CSV → Stage 4a fit on Core → Stage 4b
│                          forward apply on Extended/Holdout → write residual
│                          CSVs and decision-gate summary
└── plot_fit.py          — NEW (placeholder): parity, residual, regime-map
                            panels (§21 i–o); to be written before Stage 4
                            review.
```

The `KIN_OVERRIDE` Ref pattern is the key v6 design choice for the fit. It avoids:
1. Mutating `const` declarations in `Params` (impossible in Julia).
2. Duplicating `assembly.jl`'s `full_residual!` to thread `(j₀, α_c)` through.
3. Writing a parallel kinetics path for fitting that drifts from production.

Trade-off: the Ref is module-level state. Concurrent fits in the same Julia session would clobber each other. v6 is single-threaded by design (Newton solves are not amenable to coarse parallelism here), so this is acceptable. If parallel fits are wanted later, the override should be threaded as a function argument through `tafel_currents` instead of stored in a Ref.

### 4.2 Data

```
an_ehd/Experimental_data/
├── bloomquist_data.csv          (162 rows × 14 cols)
├── Table_S2_gap0.25mm_flow2.csv
├── Table_S3_gap0.25mm_flow6.csv
├── Table_S4_gap0.25mm_flow10.csv
├── Table_S5_gap0.5mm_flow2.csv
├── Table_S6_gap0.5mm_flow6.csv
├── Table_S7_gap0.5mm_flow10.csv
├── Table_S8_gap1.0mm_flow2.csv
├── Table_S9_gap1.0mm_flow6.csv
└── Table_S10_gap1.0mm_flow10.csv
```

Schema: `table, gap_mm, Q_total_mL_min, j_mA_cm2, phi_AN, Q_aq_mL_min, Q_org_mL_min, We_aq, We_org, FE_ADN_pct, FE_TCH_pct, FE_PN_pct, PR_ADN_kg_cm2_h, EP_ADN_kg_kWh`

Extracted from Bloomquist SI Tables S2–S10 (PNG-embedded Word tables → manual transcription → schema-normalised CSVs). 162 rows = 9 (gap × Q_total) blocks × 18 rows per block.

### 4.3 New plots (§21 i–o)

i, j: parity FE_ADN / FE_PN model vs measured, colored by gap
k, l: residual vs (j, ε_org), faceted by gap — diagnostic for bubble-physics gap on 0.25 mm
m: Lévêque δ surface over (gap, Q_total)
n: κ_eff vs ε_org, with κ_dilute and Bruggeman factor overlaid
o: We_aq–We_org regime map with Bloomquist points

---

## 5. What v6 does **not** change (preserved from v5)

- 8-species DOF layout (H⁺, OH⁻, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻, AN, ADPN, PN); Na⁺ from electroneutrality
- 9 DOFs/cell × N_mesh = 100 → 900 unknowns
- Cell-major DOF ordering, `JAC_HALFBW = 17`
- Direct Newton `(J + λI)du = −F`, λ = 10⁻¹⁰, strict L2 descent, tol 10⁻⁵
- ForwardDiff AD option with Taylor-smoothed Bernoulli for |α| < 0.01
- Convention A AN bulk concentration (per total volume)
- ε_sat = 0.0862 regime threshold; arithmetic D_mix above ε_sat
- OH⁻-pathway buffer chemistry with k₁,f = 1.4 mol/(m³·s) (Eigen 1955)
- Stage 1 / Stage 2 / Stage 3 sweeps
- All physicality checks (§13)
- Residual integrated form `J_L − J_R + S·dx`
- Bisection-based phosphate equilibrium (no Roots.jl dependency)

---

## 5b. Operational changes alongside the v6 guide

### Output directory reorganisation

`an_ehd/output/` was flat in v5 (`data/`, `plots/`, plus loose `stage*.log` files at the top). v6 reorganised it into a per-stage layout:

```
an_ehd/output/
├── cache/                        # shared solver state (217 .bin files), unchanged
├── stage1/{data,logs,plots}/
├── stage2/{data,logs,plots}/
├── stage2m/{data,logs,plots}/
├── stage4/{data,logs,plots}/     # new in v6
└── comparisons/
    └── stage2_vs_stage2m/        # was 7 stage2vs2m_*.png + summary.txt
```

This is purely an organisation change — no file content was modified, no scripts were re-run. The `comparisons/` folder holds cross-stage plots that don't belong to a single stage. Future cross-stage comparisons (e.g. Stage 4 fit vs Stage 2m baseline) go alongside.

### Sanity-table corrections

Two tables in the v6 guide were initially populated with rough hand-sketched values that disagreed with the actual implemented module computations. After smoke-testing `cell_voltage.jl` and `hydrodynamics.jl`, both tables were updated with values from the working code:

| Table | Sketched (initial v6 draft) | Computed (final v6) | Cause |
|---|---|---|---|
| §17.2 κ_dilute | 6.6 S/m | **19.1 S/m** | Sketch underweighted PO₄³⁻ (z² = 9 at pH 13.03 where most phosphate is PO₄³⁻) and Na⁺ |
| §18.2 δ_lam(0.5 mm, 2 mL/min) | 190 μm | **146 μm** | Sketch used d_h ≈ gap (slot approximation); correct is d_h = 2·gap·W/(gap+W) for a rectangular duct |

The fitting logic and the `R_series` formula are unaffected — only the order-of-magnitude reference numbers in the guide changed. The corrected κ_dilute = 19.1 S/m is *higher* than the empirical 5–10 S/m range in CONTEXT_TRANSFER §7. v6 stays with dilute theory; if Stage 4b residuals show systematic j-correlated bias, treating κ as a tier-2 fit param is the v7 escalation path.

---

## 6. Roadmap to v7 (one-line items)

1. Bubble void fraction model `ε_gas(j, gap, Q)` — Vogt 1983 or Faraday + residence-time estimate
2. `f_bubble(j, gap, Q)` enhancement on Lévêque δ — fitted single scalar or correlation
3. K_δ geometric correction on δ — promote to tier-2 fit param if needed
4. TCH species + `j_TCH ∝ c_AN³` Tafel — 9-species DOF layout (b becomes 19)
5. m_i-corrected D_eff (§4.2) re-evaluated against Bloomquist if v6 fit shows D_mix is the limiter
6. Constant-current operating mode within Newton (drop the bisection wrapper)
7. Anode Tafel breakout (η_a(j)) replacing constant V_CE if V_cell residuals demand it

---

*References for v6 additions: Newman, Electrochemical Systems 3rd ed. §11.3; Bird/Stewart/Lightfoot Transport Phenomena 2nd ed. §14.4; Lévêque, Ann. Mines 1928; Bloomquist et al. CEJ 2026 528, 172125 (and SI Tables S2–S10).*
