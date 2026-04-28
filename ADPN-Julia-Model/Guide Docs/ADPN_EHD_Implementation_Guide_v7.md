# ADPN EHD Implementation Guide — v7

**Status:** v7 written 2026-04-28. v6 remains the canonical reference for sections not touched here (governing equations, NP transport, regime-aware D_mix arithmetic, OH⁻-pathway buffer chemistry, Newton solver, continuation strategy, cell-voltage decomposition §17, hydrodynamics §18). v7 introduces the 9th species (TCH) and folds in the v6.x patches that promoted `n_ADN, n_PN` from hardcoded constants to fit parameters.

**Working directory for v7 model code:** `an_ehd_v2/` (sibling of `an_ehd/`, which is preserved as the v6.x baseline).
**Shared experimental data:** `Experimental_data/` at the project root (used by both `an_ehd/` and `an_ehd_v2/`).

> **How to read this guide.** v6 is still authoritative for sections it covers. This guide states *only what changes* — every parameter, equation, or module update v6 already documents is unchanged unless explicitly called out below. Cross-references like "v6 §5.2" point at the v6 guide.

---

## 1. Why v7 exists

The v6.x first fit (CHANGELOG_V5toV6.md §9) returned `n_ADN = 1.000 @ LB` and three of four decision gates failing — exactly the diagnostic signature for "the v6 model is missing a high-c_AN current sink." The two physically defensible interpretations were Langmuir-Hinshelwood saturation or missing TCH species, and the v6.x §9.6 Finding 1 made the call: distinguish them by adding TCH.

v7 is that step. It:

1. Adds **TCH (1,3,6-tricyanohexane)** as species 9.
2. Adds a 4th cathodic reaction (`j_TCH`) with its own kinetic params.
3. Fixes HER kinetics at the v6.x converged values (rather than literature pure-Cd values, which were 250× lower than the actual Bloomquist surface produced).
4. Folds in the v6.x patches: `n_ADN, n_PN` are fit parameters via `KIN_OVERRIDE`, and the kinetic Tafel rate uses `(c_AN/c_ref)^n_r` for arbitrary `n_r`.

Net fit dimension: **N_THETA = 9** (was 8 in v6.x, was 6 in v6).

What v7 explicitly does **not** do:
- Bubble physics (still deferred — v6 roadmap Step 6).
- V_CE / R_contact joint fit (still frozen — v6 §20.5).
- Concentrated-solution corrections to κ_dilute.
- Anode Tafel breakout.

---

## 2. New artefact: TCH species

### 2.1 Stoichiometry (per Bloomquist SI)

Confirmed from the SI: **n_e_TCH = 2** electrons per TCH molecule, with 3 AN consumed per TCH:

```
3 CH₂=CHCN + 2 H₂O + 2 e⁻ → TCH + 2 OH⁻                  (alkaline)
```

This is the same per-molecule electron count as ADPN and PN (n_e = 2), so the Faradaic-flux formulae for TCH carry the familiar `1 / (2F)` prefactor — no new convention.

Atom balance with `n_e = 2` is consistent with TCH = `C₉H₁₁N₃` (one degree of unsaturation retained from 3 AN's three C=C double bonds — likely an internal C=C or an imine; the precise structure does not affect transport modelling).

| Quantity | Symbol | Value | Notes |
|---|---|---|---|
| Electrons per TCH | n_e_TCH | **2** | Confirmed from Bloomquist SI |
| AN consumed per TCH | — | 3 | Trimer molecularity |
| OH⁻ produced per electron | — | 1 | Same convention as ADPN, PN, HER |
| Molecular formula | — | C₉H₁₁N₃ | Inferred from `n_e = 2` mass balance |
| Molecular weight | M_TCH | **161.20 g/mol** *<TODO: confirm against Bloomquist FE_TCH column derivation>* | If Bloomquist used MW=175.23 (C₉H₁₃N₃) in their FE_TCH column, residuals will carry a 175/161 = 1.087 bias — multiply through. |

### 2.2 Transport parameters

| Quantity | Symbol | Default for v7 | Source |
|---|---|---|---|
| D_TCH in aqueous phase | D_TCH,aq | **7.0×10⁻¹⁰ m²/s** *<TODO: confirm>* | Wilke-Chang scaling from ADPN's 9.0×10⁻¹⁰ (M=108) at TCH M=161: `D ∝ M⁻⁰·⁶` |
| D_TCH in organic phase | D_TCH,org | **1.2×10⁻⁹ m²/s** *<TODO: confirm>* | Typical 1.5–2× aq for the same solute in lower-viscosity organic |
| Partition coefficient | m_TCH = c_org/c_aq | **1.5** *<TODO: lab>* | TCH is more lipophilic than ADPN (m_ADPN = 0.4); estimate based on chain length / log P. **Single most uncertain transport input — highest priority to confirm.** |
| Charge | z_TCH | 0 | Neutral organic |
| Bulk concentration | c_TCH,bulk | 0 mol/m³ | Reaction product, not feedstock |

### 2.3 Kinetic parameters

| Quantity | Symbol | Initial guess | Bounds | Source |
|---|---|---|---|---|
| Standard potential | E°_TCH | **−1.30 V vs SHE** *<TODO: confirm>* | (fixed) | Same as E°_1, E°_2 (Mathison JACS 2025) unless Bui Lab has TCH-specific value |
| Exchange current density | j₀,TCH | 1.0×10⁻³ A/m² | [10⁻⁶, 10⁻¹] | Initial fit guess; matches j₀,1 / j₀,2 v5 defaults |
| Cathodic transfer coefficient | α_c,TCH | 0.5 | [0.30, 0.70] | Initial fit guess; centred between bounds |
| AN reaction order | n_TCH | 3.0 | [1.0, 3.0] | Initial fit guess matching trimer molecularity. Could relax to [0.5, 3.0] alongside n_ADN if v3 results suggest. |

---

## 3. Species and DOFs (changed from v6 §2)

### 3.1 Species ordering

```
1: H⁺           5: PO₄³⁻        9: TCH    ← NEW
2: OH⁻          6: AN
3: H₂PO₄⁻       7: ADPN
4: HPO₄²⁻       8: PN
```

TCH appended at index 9 to keep ordering stable for v6 → v7.

### 3.2 DOF layout

| | v6 | **v7** |
|---|---|---|
| n_species | 8 | **9** |
| DOFs per cell (species + φ_l) | 9 | **10** |
| Total DOFs (N_mesh = 100) | 900 | **1000** |
| Species index k → DOF index | `9·(ix−1) + k` | **`10·(ix−1) + k`** |
| φ_l index | `9·ix` | **`10·ix`** |
| Banded Jacobian halfbw | 17 | **19** |
| Sparse colors (banded FD) | 35 | **39** |

**Bandwidth derivation:** for a 1D cell-major layout with `n_dof_per_cell` DOFs per cell and a 3-point stencil (current cell + nearest neighbours), `halfbw = 2 · n_dof_per_cell − 1`. With `n_dof_per_cell = 10`, `halfbw = 19` and `n_colors = 2 · halfbw + 1 = 39`. This must be updated wherever v6 hard-coded `JAC_HALFBW = 17` or `n_colors = 35`.

---

## 4. Tafel rate expressions (rewritten v6 §5.2)

### 4.1 Four cathodic reactions

```
ADPN (R1) :  η_1 = (φ_s − φ_l) − E°_1
             j_1 = j₀,1 · (c_AN/c_ref)^n_ADN · exp(−α_c,1 · F · η_1 / RT)

PN   (R2) :  η_2 = (φ_s − φ_l) − E°_2
             j_2 = j₀,2 · (c_AN/c_ref)^n_PN  · exp(−α_c,2 · F · η_2 / RT)

HER  (R3) :  η_3 = (φ_s − φ_l) − E°_3
             j_3 = j₀,3 ·                    exp(−α_c,3 · F · η_3 / RT)

TCH  (R4) :  η_4 = (φ_s − φ_l) − E°_TCH                                        ← NEW
             j_4 = j₀,TCH · (c_AN/c_ref)^n_TCH · exp(−α_c,TCH · F · η_4 / RT)
```

Notes:
- v7 keeps E°_1 = E°_2 = −1.30 V vs SHE from v5; E°_TCH defaults to the same pending lab confirmation.
- HER (R3) has no AN concentration dependence (water reduction), unchanged from v5/v6.
- The general `(c_AN/c_ref)^n_r` form is already in `an_ehd/kinetics.jl` per v6.x; v7 just adds the 4th branch and extends the exponent tuple.

### 4.2 KIN_OVERRIDE shape (v7 extension)

The override Ref carries:

```julia
@NamedTuple{j0::NTuple{4,Float64},   # was 3-tuple in v6.x: (ADPN, PN, HER, TCH)
            ac::NTuple{4,Float64},   # was 3-tuple in v6.x: (ADPN, PN, HER, TCH)
            n ::NTuple{3,Float64}}   # was 2-tuple in v6.x: (n_ADN, n_PN, n_TCH)
```

`tafel_currents` returns a 4-tuple `(j_1, j_2, j_3, j_4)`. The HER branch ignores `n`. Default fallback (override = `nothing`) reads from Params with hard-coded `n_ADN = 2, n_PN = 1, n_TCH = 3` (so Stages 1 / 2 / 2m / 3 stay byte-identical to v5 if Params is left unchanged).

---

## 5. Faradaic flux boundary conditions (rewritten v6 §7.1)

For each reaction, current density `j_r` [A/m²] gives a product formation rate `j_r / (n_e_r · F)` [mol/m²/s] at the electrode (x = 0). Reactant fluxes follow from molecularity.

### 5.1 Updated electrode-side Neumann BCs

```
J_OH(0)   = +(j_1 + j_2 + j_3 + j_4) / F
            (every electron generates 1 OH⁻ in alkaline; same convention as v6 with TCH added)

J_AN(0)   = −(2 j_1 + j_2 + 3 j_4) / (2 F)
            (= −[ j_1/F   +  j_2/(2F)  +  3·j_4/(2F) ])
            (ADPN consumes 2 AN per molecule, PN 1, TCH 3)

J_ADPN(0) = +j_1 / (2 F)
J_PN(0)   = +j_2 / (2 F)
J_TCH(0)  = +j_4 / (2 F)                                                  ← NEW

J_other(0) = 0   (H⁺, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻ — no Faradaic source/sink at electrode)
```

The new TCH term in `J_AN(0)` is the most error-prone change. With TCH consuming 3 AN per molecule and `n_e = 2`, the AN sink from TCH is `3 · (j_4 / (2F)) = 3 j_4 / (2F)`. Verify by checking that total AN-consumption-per-electron sums correctly: `(2 e⁻ → 2 AN) + (2 e⁻ → 1 AN) + (2 e⁻ → 3 AN) = 6 e⁻ / 6 AN = 1 AN per e⁻` for ADPN, 0.5 for PN, 1.5 for TCH. This matches the integrated `J_AN(0)` formula above.

### 5.2 Bulk-side Dirichlet BCs (unchanged from v6 §7.2)

`c_TCH(δ) = 0` (no TCH in feed), `c_TCH,bulk = 0`. All other v6 §7.2 BCs unchanged.

---

## 6. Diffusivity — 9-species mixture (rewritten v6 §4)

### 6.1 Array extensions

`D_aq`, `D_org`, `m_partition`, `z_species` arrays grow from length 8 to length 9 in `an_ehd_v2/params.jl` and `an_ehd_v2/diffusivity.jl`. New entry at index 9:

```julia
D_aq[9]  = 7.0e-10   # m²/s, TCH in aqueous, <TODO: confirm>
D_org[9] = 1.2e-9    # m²/s, TCH in organic,  <TODO: confirm>
m_partition[9] = 1.5 # c_org/c_aq, TCH,        <TODO: lab>
z_species[9] = 0     # neutral
```

### 6.2 Regime-aware D_mix (no formulation changes)

The arithmetic / m_i-corrected D_mix logic from v6 §4.1, §4.2 applies unchanged to species 9. Both regimes use the same per-species `(D_aq, D_org, m)` triple — no new branches.

---

## 7. Parameter Tables — additions to v6 §9

### 7.1 §9.1 Physical constants — TCH row added

```
Species index 9: TCH (1,3,6-tricyanohexane)
   D_aq         = 7.0e-10  m²/s   <TODO: confirm>
   D_org        = 1.2e-9   m²/s   <TODO: confirm>
   m_partition  = 1.5             <TODO: lab>
   M_w          = 161.20 g/mol    <TODO: confirm vs Bloomquist FE_TCH derivation>
   z            = 0
   c_bulk       = 0   mol/m³
```

### 7.2 §9.2 Kinetic Fitting Parameters — extended

| Parameter | Initial | Bounds | Source |
|---|---|---|---|
| j₀,TCH | 1.0×10⁻³ A/m² | [10⁻⁶, 10⁻¹] | Initial fit guess |
| α_c,TCH | 0.5 | [0.30, 0.70] | Initial fit guess |
| n_TCH | 3.0 | [1.0, 3.0] | Trimer molecularity |
| E°_TCH | −1.30 V vs SHE | (fixed) | <TODO: confirm> |

### 7.3 §9.4 Partition coefficients — TCH row added

```
m_TCH = 1.5                               <TODO: lab>
```

---

## 8. Fitting Strategy (rewritten v6 §20)

### 8.1 Fit dimension: N_THETA = 9 (was 8 in v6.x, 6 in v6)

```
theta = (log10 j₀,1, log10 j₀,2, log10 j₀,TCH,
         α_c,1, α_c,2, α_c,TCH,
         n_ADN, n_PN, n_TCH)
```

**HER (`j₀,3`, `α_c,3`) is removed from the fit** and frozen at v6.x converged values:

```
j₀,3   = 2.666×10⁻⁵ A/m²
α_c,3  = 0.390
```

These are the values from `an_ehd/output/stage4v2/data/stage4a_fitted_theta.txt`. They are used because:
1. They already balance against Bloomquist's actual cell — Core HER residual was small under v6.x.
2. Pure-Cd literature values (Trasatti 1972: j₀ ≈ 10⁻⁷ A/m², α_c ≈ 0.50) are ~250× lower than what the Bloomquist surface produces, and pinning them there would inject systematic FE_HER bias into the ADPN/PN/TCH fit.
3. If v7 fit residuals show systematic FE_HER trend with j or ε_org, HER can be thawed back into the fit in v8.

### 8.2 Bounds layout

```
THETA_LB = [-6.0, -6.0, -6.0,    0.30, 0.30, 0.30,    1.0, 0.5, 1.0]
THETA_UB = [-1.0, -1.0, -1.0,    0.70, 0.70, 0.70,    3.0, 2.0, 3.0]
THETA0   = [-3.0, -3.0, -3.0,    0.50, 0.50, 0.50,    2.0, 1.0, 3.0]
                # j0    j0    j0      ac    ac    ac      n     n     n
                # ADN   PN    TCH     ADN   PN    TCH     ADN   PN    TCH
```

> **Open question for v3 design:** the v6.x fit returned `n_ADN = 1.000 @ LB`. v6.x §9.6 Finding 1 argued for relaxing this bound to 0.5 to *diagnose* whether the bound-pinning was driven by Langmuir-Hinshelwood saturation or by missing TCH. v7 *now has TCH*, so the diagnostic experiment is: keep `THETA_LB[7] = 1.0` for the first v3 fit and see whether n_ADN comes off the bound naturally. If it still pins at 1.0 with TCH present, *then* relax to 0.5 in a v4 fit and check coverage saturation.

### 8.3 Residual builder — FE_TCH added

`F` length grows from `2 · length(sel)` to `3 · length(sel)`. Per-row layout:

```
F[3n - 2] = FE_ADN_model − FE_ADN_obs   (pp)
F[3n - 1] = FE_PN_model  − FE_PN_obs    (pp)
F[3n    ] = FE_TCH_model − FE_TCH_obs   (pp)         ← NEW
```

The sigma weighting from Bloomquist's GPR (σ_ADN = 5 pp, σ_PN = 2 pp) extends to σ_TCH ≈ 3 pp *<TODO: confirm against SI>*. Only relevant if the optional weighted-loss branch is enabled.

### 8.4 Decision gates (v7 thresholds)

| Gate | v6 | v6.x | **v7 target** |
|---|---|---|---|
| Core FE_ADN RMSE | < 8 pp | < 8 pp | < 8 pp (unchanged) |
| Core FE_PN RMSE | < 5 pp | < 5 pp | < 5 pp (unchanged) |
| Core FE_TCH RMSE | — | — | **< 4 pp** (new) |
| Extended FE_ADN RMSE | < 12 pp | < 12 pp | < 12 pp |
| Holdout FE_ADN RMSE | < 15 pp | < 15 pp | < 15 pp |

The Core FE_TCH < 4 pp threshold is set tighter than FE_PN's < 5 pp because TCH FE values span a smaller range (~5–17%) than PN, so smaller absolute pp error is meaningful.

### 8.5 What if the v7 fit fails

If v7 still fails the Core ADN gate by > 5 pp, the structural issue is no longer reaction order or species count — bubble physics (Step 6 of v6 §7 roadmap) is the most likely remaining lever, with V_CE / R_contact joint refit as a smaller secondary correction.

---

## 9. Module Structure (changed from v6 §14)

All v7 files live in `an_ehd_v2/`. The structural diff from `an_ehd/`:

| File | Change |
|---|---|
| `params.jl` | `n_species = 8 → 9`; extend `D_aq`, `D_org`, `m_partition`, `z_species`, `c_bulk` arrays to length 9; add TCH constants. Bump `JAC_BLOCK = 9 → 10`, `JAC_HALFBW = 17 → 19`, `n_colors = 35 → 39`. Add `j0_TCH`, `alpha_c_TCH`, `n_TCH_default = 3.0`, `E0_TCH`, `nE_TCH = 2`. |
| `kinetics.jl` | `KIN_OVERRIDE` now carries `j0::NTuple{4,Float64}`, `ac::NTuple{4,Float64}`, `n::NTuple{3,Float64}`. `tafel_currents` returns 4-tuple `(j_1, j_2, j_3, j_4)`. New default fallback uses Params constants for all four reactions. |
| `assembly.jl` | DOF index calculations switch from `9·(ix−1)+k` to `10·(ix−1)+k` for species, `10·ix` for φ_l. Faradaic flux lines updated per §5.1. **Carefully re-check that `full_residual!` iterates `1:n_species` rather than hardcoding 8.** |
| `chemistry.jl` | `make_initial_guess(N_mesh, c_eq, eps_org)` extends to populate index 9 with `c_TCH = 0`. `bulk_concentration` returns length-9 vector. |
| `diffusivity.jl` | `D_mix` arrays become length 9; no logic changes. |
| `transport.jl` | Iterates over `1:n_species` — should already be parametric. **Audit for hardcoded 8.** |
| `solver.jl`, `mesh.jl`, `sweep_runner.jl` | No changes expected. Audit for hardcoded DOF count. |
| `fixed_j_solver.jl` | `solve_at_j` accepts `n_orders::NTuple{3,Float64}` (was 2-tuple). `j0::NTuple{4,Float64}`, `alpha_c::NTuple{4,Float64}`. Default `n_orders = (2.0, 1.0, 3.0)`. |
| `fit_kinetics.jl` | `N_THETA = 9`. `theta_to_physical` returns `(j0::NTuple{4}, ac::NTuple{4}, n::NTuple{3})` where j0[3] and ac[3] are the *frozen* HER values. `residuals!` fills `3 · length(sel)` residuals (FE_ADN, FE_PN, FE_TCH per row). Constants `J0_3_FROZEN = 2.666e-5`, `ALPHA_C3_FROZEN = 0.390` set near the top. |
| `run_stage4v3.jl` | NEW. Mirrors `run_stage4v2.jl` but writes to `an_ehd_v2/output/stage4v3/`. Saves fitted j₀,TCH, α_c,TCH, n_TCH alongside the others. |
| `run_stage4.jl`, `run_stage4v2.jl` | Both gated with `error()` to prevent overwriting v6 / v6.x baselines. |
| `analyze_stage4.jl` | Update to read 9-key fitted-theta and emit `stage4v3_diagnostic.csv` with FE_TCH column added. |
| `plot_stage4_*.py` | Updates: target `an_ehd_v2/output/stage4v3/`, add FE_TCH parity panel, FE_HER residual now uses corrected `100 − FE_ADN − FE_TCH − FE_PN` denominator (TCH no longer absorbed silently). |

### 9.1 Path note — Experimental_data moved

`Experimental_data/` was moved to the project root for shared use across model versions. v7 code paths must update from:

```julia
joinpath(@__DIR__, "Experimental_data", "bloomquist_data.csv")  # v6 (relative to an_ehd/)
```

to:

```julia
joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")   # v7
```

Affected files: `run_stage4v3.jl`, `analyze_stage4.jl`, any plot script that reads `bloomquist_data.csv`. Audit before running.

---

## 10. Re-run requirements

The 9-species DOF layout invalidates every v6 cache. Stages 1, 2, 2m, 3 must be re-run before Stage 4v3. Order:

1. **Stage 1** (single-point baseline) — verify the v7 model converges at all. Fast.
2. **Stage 2** (arithmetic ε_org sweep) — confirm trends match v6.
3. **Stage 2m** (m_i-corrected D_eff sweep) — same. Optional unless m_i upgrade is being re-evaluated.
4. **Stage 3** (warm-start cache builder) — required before Stage 4v3 to keep the LM fit time tractable. Builds V = −1.0 V states for every unique `(gap, Q_total, ε_org)` key in the Bloomquist data.
5. **Stage 4v3** — the new fit.

Total wake compute estimate: ~30 min for Stage 1+2+2m, ~15 min Stage 3, ~75 min Stage 4v3 = **~2 hours wake compute**. Use `caffeinate -dimsu -w <julia_pid>` and stay on AC; lid-closure on battery will pause the process even with caffeinate (the `-s` flag is documented as AC-only).

---

## 11. Plots — required updates

In addition to the v6 plot panels (i)–(o), v7 needs:

- **`plot_stage4v3_parity.py`** — adds 4-panel composite: FE_ADN, FE_PN, **FE_TCH**, V_cell parity (4 panels, 2×2).
- **`plot_stage4v3_residuals.py`** — adds FE_TCH residual vs (j, ε_org), faceted by gap.
- **`plot_stage4v3_3d_surfaces.py`** — adds 4th panel (FE_TCH) to the regime map.

Old plot scripts for stage4 and stage4v2 should be left intact (they target frozen baselines).

---

## 12. v7-specific pitfalls

1. **Bandwidth bug.** Anywhere v6 hardcoded `JAC_HALFBW = 17` (e.g. assembly.jl, solver.jl) must update to 19. Banded FD jacobian colors must update to 39.

2. **`n_species` hardcoded as 8.** Search the codebase for literal `8` in array indexing — there should be very few. The v6 code is mostly parametric over `Params.n_species`, but verify every loop bound.

3. **DOF index off-by-cell.** The `9·(ix−1) + k` → `10·(ix−1) + k` change is mechanical but easy to miss in one obscure place. After each module is updated, run `test_smoke.jl` to verify Newton converges from a cold IC.

4. **AD type stability through `tafel_currents` returning 4 instead of 3 values.** ForwardDiff specialises on tuple length. After updating `assembly.jl`'s caller, rebuild types from a fresh Julia session.

5. **FE_TCH residual bias from MW mismatch.** If Bloomquist computed FE_TCH using MW = 175.23 (C₉H₁₃N₃) but v7 uses MW = 161.20 (C₉H₁₁N₃), the model's FE_TCH will be biased low by 161/175 ≈ 0.92. This is a 1.087× systematic correction *<TODO: confirm by reading the SI's FE_TCH column definition>*.

6. **HER frozen-value mismatch with Params.** v6 Params has `j0_3 = 1e-5`, but v7 freezes at v6.x's fitted `j0_3 = 2.666e-5`. The discrepancy must be either (a) updated in v7 Params, or (b) applied via a permanent KIN_OVERRIDE in the fit driver. Pick one — mixing causes confusion.

---

## 13. v7 → v8 roadmap

After v7 lands and produces a Stage 4v3 fit, the remaining structural physics is bubbles. v8 priorities, ordered:

1. **`f_bubble(j, gap, Q)` enhancement on δ_lam** (v6 §7 Step 6c) — closes the Holdout (0.25 mm gap) gap.
2. **Bruggeman void factor on κ_eff** (v6 §7 Step 6a) — secondary effect; helps V_cell parity.
3. **`ε_gas(j, gap, Q)` model** (v6 §7 Step 6b) — Faraday + residence-time as the first-cut, scalar `c_bubble` fit param if needed.
4. **(conditional) `n_ADN` lower bound relaxed to 0.5** — only if v7 still pins n_ADN at LB after TCH lands. That would point at coverage saturation, not the missing-species hypothesis.
5. **(conditional) Re-fit V_CE, R_contact** — only after bubbles are in.

---

## 14. What v7 does NOT change (preserved from v6)

Everything in v6 not enumerated above. In particular:

- Governing equations (§3), Scharfetter-Gummel flux (§3.3), buffer chemistry (§6), bulk-side BCs (§7.2), three-parameter sweep philosophy (§8), Newton solver and continuation (§10.3–§10.7), residual integrated form (§10.7), solution caching scheme (§11), physicality checks (§13), potential referencing (§16), cell-voltage decomposition (§17), hydrodynamics + Lévêque (§18).
- The Core / Extended / Holdout row selection (v6 §20.1) and the `j ≤ 190 mA cm⁻²`, `ε_org ≥ 0.04`, `gap ∈ {0.5, 1.0} mm` filters.
- V_CE = 1.7 V and R_contact = 1×10⁻⁴ Ω·m² defaults remain frozen for v7.

---

*References: Bloomquist et al. CEJ 2026, 528, 172125 (and SI Tables S2–S10) for the experimental data and FE column definitions; Mathison JACS 2025 for E°_1 / E°_2; Trasatti, J. Electroanal. Chem. 1972 for HER literature comparison; v6 implementation guide (`ADPN_EHD_Implementation_Guide_v6.md`) for unchanged sections; v5→v6 changelog (`ADPN_EHD_Implementation_Guide_CHANGELOG_V5toV6.md`) §9 for v6.x fit results.*
