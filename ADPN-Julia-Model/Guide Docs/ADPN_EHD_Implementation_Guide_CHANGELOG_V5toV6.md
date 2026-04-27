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

## 7. Concrete next-steps plan (post-first-fit)

The first end-to-end Stage 4 fit (run on 2026-04-27 with default kinetic guesses) plateaued near combined RMSE ≈ 10 pp on the Core subset. That's the floor of what the v6 model can represent without additional physics. The next steps below are ordered by effort × value, each scoped enough to act on without re-deriving the design.

### Step 1 — Diagnostic plots from the current fit (today, ~1 hour)

No model re-runs needed; everything reads from the residual CSVs and `analyze_stage4.jl` output.

**`plot_stage4_parity.py`** — four-panel matplotlib figure:
- (a) FE_ADN parity model vs Bloomquist, 162 points colored by gap, slope-1 dashed reference, RMSE annotation per subset.
- (b) Same for FE_PN.
- (c) Residual (model − obs) vs j, faceted by gap. Random scatter on Core ⇒ kinetics-fit-OK; systematic positive bias on Holdout ⇒ bubble physics matters.
- (d) Residual vs ε_org, faceted by gap. Systematic vs ε_org ⇒ D_mix arithmetic-mean is too weak and m_i correction (§4.2) is needed.

**`plot_stage4_3d_surfaces.py`** — recreates Bloomquist Fig. 5:
- Three panels (one per gap: 0.25, 0.5, 1.0 mm).
- Axes: log10(We_aq), log10(We_org), j; color = FE_ADN.
- Overlay model predictions on top of experimental data.
- The shape test: does the model put the FE_ADN > 70% region at high We_org / low We_aq like the paper, *independent of magnitude*? If yes, the kinetics-form is right and we just need to close the magnitude gap. If no, structural model error remains.

**`plot_stage4_v_cell_parity.py`** — V_cell_pred vs V_cell_obs (back-derived from EP/PR/j as `V_cell = PR_ADN / (EP_ADN · j_A_cm2)`). Tests whether frozen V_CE = 1.7 V and R_contact = 1×10⁻⁴ Ω·m² are reasonable. Slope-1 with small offset ⇒ fine. Slope ≠ 1 or large bias ⇒ V_CE / R_contact need fitting in v7.

### Step 2 — Loosen `tol_rel` and build Stage 3 cache (this week, ~1 hour total)

**Tol loosen:** change `tol_rel` default in `lm_fit` from 1e−4 to 1e−2. The first fit's iter 9 would have terminated at 1.31% drop; this saves ~10 LM iterations and ~30 min wall time per fit cycle. Keep 1e−4 as an option for "final" fits once the model physics stabilises in v7.

**`run_stage3.jl`** (≈80 lines):
- Loop the unique `(gap, Q_total, ε_org)` keys across Core ∪ Extended ∪ Holdout (~90 keys).
- For each key: derive δ from Lévêque, build mesh, run v5 `bootstrap!` (α_buf 0→1 ramp + α_kin 1e−6→1 ramp at V = −1.0 V).
- Write converged DOF vector to `output/cache/s_eo<ε>_d<δμm>_V-1.000.bin` (matches v5 binary cache filename convention).
- Also emit `output/stage3/data/stage3_warmstart_index.csv` listing `(key → cache filename)`.

Add file-loader path to `FitContext.build_context` so every Stage 4 invocation auto-reads the cache. After Stage 3 build, fresh `run_stage4.jl` re-runs skip the cold V-walk entirely. Cost-payback after ~3 re-runs.

### Step 3 — Reaction order n₁, n₂ as fit params (next, ~30 min)

Currently `j_1 ∝ c_AN²` (ADPN) and `j_2 ∝ c_AN¹` (PN) are hard-coded in `kinetics.jl`. Promote to fit params:

```julia
j_1 = j0_1 · (c_AN/c_ref)^n_1 · exp(-α_c1 · F · η1 / RT)
j_2 = j0_2 · (c_AN/c_ref)^n_2 · exp(-α_c2 · F · η2 / RT)
```

| Effort | Where to edit |
|---|---|
| ~30 min | Add `N1, N2 ∈ Ref` alongside `KIN_OVERRIDE` in `kinetics.jl`; modify `tafel_currents`. Add `n_1, n_2` to `theta` vector and `THETA_LB/UB/0` in `fit_kinetics.jl`. |

Bounds: physically `n_1 ∈ [1, 3]` (Tafel dimerization mechanism debated for 60 years), `n_2 ∈ [0, 2]`. Initial guesses `(n_1, n_2) = (2.0, 1.0)` match v6's hard-coded defaults.

Fit dimension: 6 → 8. Still very overdetermined (96 Core residuals). **Watch for** the optimiser drifting `n_1` toward 0.5 to fit FE shape vs ε_org — that's a sign D_mix is wrong, not that ADPN kinetics is half-order.

### Step 4 — Experimental error weights (~20 min)

Bloomquist's GPR surrogate uses noise σ_FE,ADN = 5 pp, σ_FE,PN = 2 pp (SI §"Building Surrogate Models"). v6 weights every residual equally. Switch to:

```
weighted_loss = Σ ((FE_model_i - FE_obs_i) / σ_i)²
```

| Effort | Where |
|---|---|
| ~20 min | Add `sigma_FE_ADN_pp::Float64` and `sigma_FE_PN_pp::Float64` to `BloomquistRow` (or hard-code in `residuals!`). Bloomquist doesn't tabulate per-row σ, so use the SI's global values. |

Doesn't change which rows pass the §20.4 gates, but rebalances the fit toward giving noisy rows less influence. Worth doing if reviewers ask; not critical for the test-drive.

### Step 5 — TCH species (~3 hours code + ~1 hour fit)

TCH (1,3,6-tricyanohexane) is currently absent from v6's 8-species model but is 5–17% of total FE in Bloomquist's data — that current is being absorbed into the model's FE_HER residual. Adding it:

| What | Effort |
|---|---|
| Update `Params.n_species`, `D_aq`, `D_org`, `m_partition`, `z_species` to 9 species | ~10 min |
| Update `chemistry.jl::make_initial_guess` and `bulk_concentration` | ~10 min |
| Add `j_TCH = j0_TCH · (c_AN/c_ref)^3 · exp(-α_c,TCH · F · η / RT)` to `kinetics.jl`. Stoichiometry: 3 AN + 6 e⁻ → TCH + 6 OH⁻, so n_e = 6 (or rescale Tafel by n_e/2 to keep convention). | ~15 min |
| Update `assembly.jl` Faradaic flux BCs: `J_AN(0) -= 3·j_TCH/(n_e·F)`, `J_TCH(0) = +j_TCH/(n_e·F)`, `J_OH(0)` already covers `+(j_total)/F` because TCH produces OH⁻ at the same rate per electron | ~15 min |
| Bump `JAC_BLOCK = 9 → 10`, `JAC_HALFBW = 17 → 19`, `n_colors = 39` | ~5 min |
| Add `FE_TCH` residual term to `residuals!`; fit dim grows by 2 (j₀_TCH, α_c,TCH) | ~20 min |
| Re-run Stage 4 | ~1 hour |

Easy in absolute terms, just touches several files. Should produce meaningful improvement because TCH FE is currently mismatched in the v6 fit (the model puts that current into HER instead).

### Step 6 — Bubble physics (~6–10 hours, the biggest lift)

The headline Bloomquist finding ("bubble-induced convection dominates") is exactly what v6 deferred. Three sub-steps, do in this order:

**Step 6a — Bruggeman void factor on κ_eff** (~10 min):
```
κ_eff = κ_dilute · (1 - ε_org)^1.5 · (1 - ε_gas)^1.5
```
in `cell_voltage.jl`. Free physics. No fit params (assuming we have an `ε_gas` model).

**Step 6b — `ε_gas(j, gap, Q)` model** (~3 hours):
Two paths:

| Path | Pros | Cons |
|---|---|---|
| **Faraday + residence-time** (physics-based) | Zero new fit params; predictive | Requires bubble-detachment radius and gas hold-up correlations — both empirical; brittle |
| **Single fit param `c_bubble`** (correlation `f = 1 + c·j^β`) | Robust; one new param | Less interpretable; absorbs other missing physics |

Recommended: start with Faraday `Q_H₂ = j·A/(2F)` and a fixed `τ_residence = L_channel/v_super`, giving `ε_gas = Q_H₂·τ/V_channel` × correction factor. If 6a+6b alone close the 0.25 mm holdout RMSE gap, stop. If not, promote to fit param.

**Step 6c — `f_bubble(j, gap, Q)` enhancement on δ_lam** (~30 min):
```
δ_actual = δ_lam · f_bubble(j, gap, Q)
```
in `hydrodynamics.jl`. Vogt 1983 or fit-params correlation. The stub is already in `hydrodynamics.jl` as a commented-out `delta_actual` function.

**Total v7 effort**: ~6–10 hours assuming empirical correlation form. Doubles if all three go physics-based.

### Step 7 — Stage 4c: joint refinement on (V_CE, R_contact) (~1 hour, defer until v7)

Once bubble physics is in (Step 6), revisit `(V_CE, R_contact)`. Need to think carefully about loss weighting because V residuals (in V) and FE residuals (in pp) are different units; weight by `(1/σ²)` to make them commensurate. Defer until then because:
- The frozen V_CE = 1.7 V might already be ~within 0.2 V of optimal — fitting it only buys tiny FE improvement.
- Bubble physics affects R_series too — V_CE/R_contact you'd fit in v6 would be wrong in v7 once bubbles land. Better to wait for the post-bubble landscape.

### Stop condition

**Don't do all of these.** Stop adding parameters when:
1. Residuals look random vs (j, ε_org, gap) on the parity plots.
2. Combined RMSE on Core is in the noise floor (~σ_FE,ADN ≈ 5 pp).
3. Adding the next param doesn't reduce parity-plot scatter visibly.

Beyond that, you're fitting noise. Stop and start interpreting the fit physically.

### Recommended sequence in dependency order

```
Step 1 (plots) ─┐
                ├─→ Step 3 (reaction orders) ─┐
Step 2 (Stage 3 + tol) ─┘                     ├─→ Step 5 (TCH) ─→ Step 6 (bubbles) ─→ Step 7 (V_CE/R_contact)
                                              │
Step 4 (weights) ─────────────────────────────┘  (optional, anytime after Step 1)
```

Steps 1–4 are v6.x patches. Step 5 onward is the v7 boundary.

---

## 8. First-fit results (2026-04-27)

The first end-to-end Stage 4 run completed on 2026-04-27 with the loosened LM gates (`tol_rel = 1e-2`, `lambda_stuck = 1e3`) and the Stage 3 warm-start cache (126 / 126 unique transport keys built cleanly). The intent was a test-drive of the fitting pipeline, not a publishable fit. Results below are diagnostic of v6's structural model limits, not final fitted parameters.

### 8.1 Run summary

| | Value |
|---|---|
| Core rows fit | 48 (gap ∈ {0.5, 1.0} mm, j ≤ 190 mA/cm², ε_org ≥ 0.04) |
| Extended rows (forward-applied, no re-fit) | 42 (gap ∈ {0.5, 1.0} mm, ε_org ≥ 0.04, any j) — Extended-only = Extended ∖ Core |
| Holdout rows (forward-applied, no re-fit) | 43 (gap = 0.25 mm, ε_org ≥ 0.04) |
| LM iterations to converge | 11 (rel-loss-drop < 1e-2 gate) |
| Wall time (with Stage 3 cache + loose gates) | ≈ 12 min |
| All decision-gate checks | 1 / 4 PASS |

### 8.2 Fitted kinetic parameters

| Param | Initial | Bounds | Fitted | Pinned at bound? |
|---|---|---|---|---|
| **j₀,1 (ADPN)** | 1.0×10⁻³ | [10⁻⁶, 10⁻¹] A/m² | **1.000×10⁻⁶** | ⚠️ at lower bound |
| j₀,2 (PN) | 1.0×10⁻³ | [10⁻⁶, 10⁻¹] A/m² | 6.166×10⁻³ | free |
| j₀,3 (HER) | 1.0×10⁻⁵ | [10⁻⁸, 10⁻³] A/m² | 9.533×10⁻⁵ | free |
| **α_c,1 (ADPN)** | 0.5 | [0.3, 0.7] | **0.700** | ⚠️ at upper bound |
| α_c,2 (PN) | 0.5 | [0.3, 0.7] | 0.332 | near lower, not pinned |
| α_c,3 (HER) | 0.4 | [0.3, 0.5] | 0.308 | near lower, not pinned |

### 8.3 Decision-gate scoreboard (§20.4)

| Gate | Result | Threshold |
|---|---|---|
| Core FE_ADN RMSE | **13.50 pp FAIL** | < 8 pp |
| Core FE_PN RMSE | **4.84 pp PASS** | < 5 pp |
| Extended FE_ADN RMSE | **14.25 pp FAIL** | < 12 pp |
| Holdout FE_ADN RMSE | **34.85 pp FAIL** | < 15 pp |

Three of four gates failed exactly along the modes anticipated in §20.4 and §6 (v7 roadmap). The PN gate passed cleanly.

### 8.4 V_cell parity (post-fit, `analyze_stage4.jl`)

V_cell back-derived per row from `V_cell_obs = 1000 · PR_ADN / (EP_ADN · j_A_cm2)` and compared with `V_cell_pred = V_CE + |V_cathode_SHE| + j · R_series` using frozen `V_CE = 1.7 V`, `R_contact = 1×10⁻⁴ Ω·m²`:

| Subset | n | median(V_pred − V_obs) | MAE | bias |
|---|---|---|---|---|
| Core | 48 | +0.275 V | 0.376 V | +0.233 V |
| Extended-only | 42 | −0.405 V | 0.445 V | −0.405 V |
| Holdout (0.25 mm) | 43 | +0.353 V | 0.436 V | +0.230 V |
| Excluded (ε<0.04) | 20 | −0.149 V | 0.741 V | −0.090 V |

> **Caveat on V_cell_obs.** The 1000× factor (kW→W) was missed in the first analyze_stage4 run; numbers above are from the corrected re-run. PR_ADN, EP_ADN, j each carry independent measurement noise → V_cell_obs has compounded uncertainty estimated at 0.2–0.4 V.

### 8.5 Interpretation — three findings

**Finding 1 — ADPN kinetics form is structurally too rigid.**
LM pushed both j₀,1 (to lower bound) and α_c,1 (to upper bound) simultaneously, yet still left FE_ADN RMSE = 13.5 pp on Core. Bound-pinned parameters in opposite directions are the classical signature of a model that can't represent the data shape. The fixed `c²` AN dependence in `j_1 = j₀,1 · (c_AN/c_ref)² · exp(...)` is the suspect — Mathison (JACS 2025) and the broader 60-year ADPN literature show evidence the effective reaction order can sit anywhere in [1, 3] depending on potential and surface coverage. **Action: promote `n_1` to a fit param** (Step 3 of §7 roadmap).

**Finding 2 — PN kinetics fits cleanly. Don't touch.**
FE_PN Core RMSE = 4.84 pp passes the < 5 pp gate. j₀,2 free, α_c,2 not bound-pinned. The `c¹` dependence for PN appears adequate. No edits needed in v7.

**Finding 3 — Bubble physics is the dominant missing term, but only on transport (FE), not on V_cell.**
The Holdout − Extended-only delta on FE_ADN is +20.6 pp (34.85 − 14.25), with only the gap differing (0.25 mm vs 0.5 / 1.0 mm). However, the Holdout V_cell bias (+0.23 V) is *similar* to the Core bias, not catastrophic. **Diagnostic implication:** the v6 missing-physics is dominated by **bubble-induced convection** (which enhances mass transfer at small gaps, fixing FE) and **not** by bubble void blocking of κ_eff (which would also push V_cell up). This argues that `f_bubble(j, gap, Q)` enhancement on δ_lam (Step 6c of §7 roadmap) is more important than the Bruggeman correction on κ_eff (Step 6a). Both should land in v7, but if priorities force a choice, do 6c first.

**Bonus finding — V_CE / R_contact defaults look fine.**
MAE on V_cell ≈ 0.4 V across all subsets is within compounded measurement noise of V_cell_obs. The Core (+0.28 V over-prediction) vs Extended-only (−0.41 V under-prediction) sign flip is consistent with R_contact being slightly small at high j, but the magnitude is too small to fit reliably against the back-derived V_cell. Keep frozen in v7 until the bubble model lands.

### 8.6 Output artefacts (commit-worthy diagnostic CSVs)

```
an_ehd/output/stage4/data/
├── stage4a_fitted_theta.txt          (354 B)
├── stage4a_core_residuals.csv        (3.0 KB, 48 rows)
├── stage4b_extended_residuals.csv    (5.7 KB, 90 rows)
├── stage4b_holdout_residuals.csv     (2.9 KB, 45 rows)
└── stage4_diagnostic.csv             (21 KB, 162 rows × 21 cols — full enriched output incl. V_cell)
```

### 8.7 Confirmed v7 priority order

The first fit's failure pattern unambiguously orders the v7 work that §7 listed without empirical justification:

1. **Reaction order n_1 as fit param** (Step 3) — most direct fix for the bound-pinned ADPN kinetics. Cheap (~30 min code).
2. **Bubble convection on δ_lam** (Step 6c) — biggest residual remaining after Step 1, dominates the 0.25 mm gap holdout.
3. **Bubble Bruggeman on κ_eff** (Step 6a) — secondary effect; helps V_cell parity but not the FE residuals.
4. **TCH species** (Step 5) — doesn't show up clearly in current residuals because TCH's ~10% FE fraction is being absorbed into model FE_HER without distorting FE_ADN/FE_PN much. Lower priority than originally suggested.
5. (still deferred) V_CE / R_contact joint fit (Step 7) — only meaningful after bubbles land.

---

*References for v6 additions: Newman, Electrochemical Systems 3rd ed. §11.3; Bird/Stewart/Lightfoot Transport Phenomena 2nd ed. §14.4; Lévêque, Ann. Mines 1928; Bloomquist et al. CEJ 2026 528, 172125 (and SI Tables S2–S10).*
