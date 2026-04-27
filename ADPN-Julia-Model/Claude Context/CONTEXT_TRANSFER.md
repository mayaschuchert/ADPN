# ADPN EHD Model — Context Transfer Document

**Project:** 1D planar Nernst-Planck model for acrylonitrile (AN) electrohydrodimerization (EHD) to adiponitrile (ADPN) on a Cd cathode.
**Lab:** Bui Lab, NYU Tandon
**Working directory:** `C:\Users\BuiLab\Documents\ADPN-Julia-Model\`
**All model code lives in:** `an_ehd/`

---

## 1. What Has Been Built and Is Working

Stages 1, 2, and 2m are complete — all CSV output exists.

| Stage | Description | Output tag | ε_org values |
|-------|-------------|------------|--------------|
| 1 | Single point, baseline | `stage1` | 0.02 |
| 2 | Arithmetic D_mix sweep | `stage2` | 0.02, 0.05, 0.09, 0.15, 0.25, 0.30 |
| 2m | m_i-corrected D_eff | `stage2m` | same |

Records CSVs: `an_ehd/output/data/<tag>_records_eo<eps>_d50.csv`
Profile CSVs: `an_ehd/output/data/<tag>_profile_eo<eps>_d50_V<V>.csv`
Comparison plots: `an_ehd/output/plots/stage2vs2m_*.png`

---

## 2. File Structure and What Each File Does

```
an_ehd/
├── ADPN_EHD.jl         — master module; includes all submodules
├── params.jl           — all physical constants
├── mesh.jl             — geometric grading mesh
├── diffusivity.jl      — regime-aware D_mix (arithmetic / m_i-corrected switch)
├── chemistry.jl        — phosphate equilibrium, buffer sources, bulk conc, initial guess
├── kinetics.jl         — Tafel currents j1, j2, j3
├── transport.jl        — Scharfetter-Gummel flux
├── assembly.jl         — full_residual!, DOF layout
├── solver.jl           — Newton solver + V-continuation + log-j continuation
├── sweep_runner.jl     — run_one_sweep pipeline (mesh → equilibrium → bootstrap → V-cont → export)
├── run_stage1.jl       — thin wrapper: single sweep ε=0.02, δ=50 μm
├── run_stage2.jl       — ε_org loop sweep, arithmetic D_mix
├── run_stage2m.jl      — same loop, m_i-corrected D_eff
├── plot_results.py     — profile + polarization plots for stage1
├── plot_stage2.py      — multi-ε_org overlay plots
└── plot_stage2_comparison.py  — arithmetic vs m-corrected side-by-side
```

---

## 3. Key Physical Model Choices

### Species ordering (8 species)
1:H⁺ 2:OH⁻ 3:H₂PO₄⁻ 4:HPO₄²⁻ 5:PO₄³⁻ 6:AN 7:ADPN 8:PN

### DOF layout — cell-major, 9 DOFs per cell
- species k ∈ 1..8 → DOF index = 9·(ix−1) + k (stored as log concentration)
- φ_l → DOF index = 9·ix
- Total DOFs: 9·N_mesh = 9·100 = 900

### Boundary conditions
- x = 0 (electrode): Faradaic flux BCs for species; V is the solid-phase potential (scan variable)
- x = δ (bulk): Dirichlet on log(c_bulk) and φ_l = 0

### Tafel kinetics (`kinetics.jl`)
Three reactions — ADPN (R1), PN (R2), HER (R3):
```
η_r = (V − φ_l_surface) − E0_r
j1 = j0_1 · (c_AN/c_ref)² · exp(−α_c1 · F · η1 / RT)   [ADPN, 2nd order in AN]
j2 = j0_2 · (c_AN/c_ref)  · exp(−α_c2 · F · η2 / RT)   [PN,   1st order]
j3 = j0_3 ·                 exp(−α_c3 · F · η3 / RT)   [HER,  no AN dep.]
```
Parameters in `params.jl`:
- E0_1 = E0_2 = −1.30 V vs SHE, E0_3 = −0.83 V vs SHE
- α_c1 = 0.5, α_c2 = 0.5, α_c3 = 0.4
- j0_1 = j0_2 = 1e-3 A/m², j0_3 = 1e-5 A/m²
- c_ref = 1000 mol/m³

### Electrode BC fluxes
```
J_OH(0)   = +(j1 + j2 + j3) / F       # OH⁻ produced
J_AN(0)   = −(2j1 + j2) / (2F)        # AN consumed
J_ADPN(0) = +j1 / (2F)                # ADPN produced
J_PN(0)   = +j2 / (2F)                # PN produced
```

### Scharfetter-Gummel flux (`transport.jl`)
```
α = z_i · F · (φ_R − φ_L) / (RT)
B(α) = α / (e^α − 1)      — Bernoulli function
N_i = −D/dx · [B(α)·c_R − B(−α)·c_L]
```
**Critical:** Taylor expansion for |α| < 0.01 (not a hard branch) to keep ForwardDiff AD derivatives smooth through α = 0.

### Buffer chemistry (`chemistry.jl`) — OH⁻-pathway
```
R1: H₂O ⇌ H⁺ + OH⁻          k1f=1.4 mol/(m³·s), k1r=1.4×10⁸ m³/(mol·s)
R2: OH⁻ + H₂PO₄⁻ ⇌ HPO₄²⁻  k2f=1×10⁵ m³/(mol·s), k2r=15.87 s⁻¹
R3: OH⁻ + HPO₄²⁻ ⇌ PO₄³⁻   k3f=2×10³ m³/(mol·s), k3r=44444 s⁻¹
```
Sources: S_H = r1; S_OH = r1−r2−r3; S_H2PO4 = −r2; S_HPO4 = r2−r3; S_PO4 = r3

### AN bulk concentration — Convention A
```
ε_org < ε_sat (≈0.0862): c_AN = ε_org · (ρ_AN/M_AN)    [single-phase, linear]
ε_org ≥ ε_sat:            c_AN = C_AN_SAT ≈ 1310 mol/m³  [two-phase, saturated]
```
`ε_sat = C_AN_SAT / MOLAR_DENSITY_AN = 1310/15191 ≈ 0.0862`

### Diffusivity (`diffusivity.jl`)
```
ε_org < ε_sat:   D_mix = D_aq[i]                        (single-phase)
ε_org ≥ ε_sat, arithmetic:   D_mix = ε·D_org + (1−ε)·D_aq
ε_org ≥ ε_sat, m_corrected:  D_mix = ε·m_i·D_org + (1−ε)·D_aq
```
Switch: `set_D_formulation!(:arithmetic)` or `set_D_formulation!(:m_corrected)`

Partition coefficients m_i: AN=11.59, ADPN=7.72, PN=11.59, ions=0

Aqueous diffusivities (m²/s): H⁺ 9.31e-9, OH⁻ 5.27e-9, H₂PO₄⁻ 0.846e-9,
HPO₄²⁻ 0.690e-9, PO₄³⁻ 0.610e-9, AN 2.30e-9, ADPN 1.50e-9, PN 2.30e-9

Organic diffusivities (m²/s): AN 6.00e-9, ADPN 3.90e-9, PN 6.00e-9, ions 0

---

## 4. Solver Details

### Newton solver (`solver.jl: newton_solve!`)
- Direct `(J + λI)du = −F` with **fixed λ = 1×10⁻¹⁰**
- Do NOT use adaptive LM / normal equations — tested, picks wrong branch
- Jacobian half-bandwidth **b = 17** (cell-major 9×9 blocks: b = 2·9−1 = 17)
- FD uses 2·17+1 = 35 color groups
- ForwardDiff AD option: `jacobian_mode = :ad` (use for high-V / ill-conditioned regions)
- Per-DOF step clamping: Δu_log ≤ 5.0, Δφ_l ≤ 0.015 V
- Backtracking: strict L2 descent (no slack)
- tol = 1×10⁻⁵ (|F|∞)

### Continuation strategy (`solver.jl: newton_continuation`)
- V sweeps −1.0 → −2.5 V (Phase A: AD Jacobian, adaptive ds_init=0.05)
- If Phase A stalls before −2.5 V: Phase B log-j continuation from V≈−2.0 restart
- AD Jacobian is the only option that reliably reaches −2.5 V (FD fails ~−1.8 V)

### Bootstrap (`sweep_runner.jl: bootstrap!`)
- α_buf ramp: 0 → 1 in 10 uniform steps (α_kin=0 throughout)
- α_kin ramp: 1×10⁻⁶ → 1 geometric ×2, 21 steps (α_buf=1)
- Ensures smooth entry into fully physical solution

---

## 5. Key Bugs Fixed (Do Not Reintroduce)

| Bug | Symptom | Fix |
|-----|---------|-----|
| Roots.jl blocked by Windows WDAC | import error | Replaced with inline `_bisect()` in chemistry.jl |
| Faraday constant name collision | `F` shadowed `Params.F` in assembly | Renamed residual arg from `F` to `res` |
| Wrong Jacobian half-bandwidth | Singular J, no convergence | JAC_HALFBW = 17 (not 9) |
| Adaptive LM picked wrong branch | Unphysical steady state | Reverted to fixed λ=1e-10 |
| Merit slack 1.1 stalled | Progress accepted at noise level | Strict descent (no slack) |
| Hard SG branch at α=0 | ForwardDiff gives dJ/dφ=0 at equilibrium | Taylor expansion for |α|<0.01 |
| c_AN_bulk wrong formula | 1379 mol/m³ at ε=0.05 (above sat) | Convention A: ε·ρ/M |
| D_mix not regime-aware | Arithmetic mean even below ε_sat | ε_sat threshold in D_mix |
| Stage 1 at ε=0 | log(0) = −Inf | Changed to ε_org = 0.02 |
| Buffer H⁺-pathway vs OH⁻-pathway | Wrong kinetic timescales | OH⁻-pathway k values, S_OH = r1−r2−r3 |

---

## 6. Current Model Results Summary

At δ = 50 μm (Stage 2, arithmetic D_mix):
- Peak FE_ADPN ≈ 35–45% across ε_org sweep
- Bloomquist experimental target: 73–80% FE at j = 70–300 mA/cm²

At δ = 50 μm (Stage 2m, m_i-corrected D_eff):
- Peak FE_ADPN ≈ 38–49% (≈ +3–5 pp improvement over arithmetic)
- Still well below Bloomquist 73–80%

**Gap to experiment:** The model underpredicts FE by ~30 pp. The dominant lever is kinetics (j0, α_c) rather than transport. The m_i correction is too small to explain the discrepancy.

Voltage range: Model solved at V = −1.0 to −2.5 V vs SHE.
Experimental Bloomquist j = 70–300 mA/cm² corresponds to V_cell ≈ 2.5–3.5 V (full cell), which maps to V_cathode ≈ −1.2 to −2.2 V after subtracting anode potential + ohmic drop.

---

## 7. Next Step the User Wants to Do (NOT YET STARTED)

The user wants to incorporate the Bloomquist et al. (CEJ 2026) experimental dataset and fit the model to it. The framework discussed (but not yet implemented) is:

### Cell-voltage decomposition
```
V_cathode_SHE = V_CE − V_cell_mag + j · R_series
R_series = gap / σ_e + R_contact   [Ω·m²]
```
where:
- V_cell_mag: full-cell voltage magnitude (positive, from Bloomquist paper)
- V_CE: effective anode potential vs SHE (OER on SS + overpotential), fit range 1.5–2.0 V
- σ_e: electrolyte conductivity [S/m], fit range 5–10 S/m
- R_contact: contact resistance [Ω·m²], fit range 0.5×10⁻⁴–2×10⁻⁴ Ω·m²

### Lévêque boundary layer
```
d_h = 2·gap·W / (gap + W)   [W = channel width, 4 mm]
Re = v·d_h / ν,   v = Q / (gap·W),   ν = 1×10⁻⁶ m²/s
Sh = 1.85·(Re·Sc·d_h/L)^(1/3)   [Lévêque, laminar developing BL]
δ_lam = d_h / Sh
δ_actual = K_δ · δ_lam          [K_δ = geometry correction, fit range 0.3–3.0]
```

### Weber number
```
We = ρ_aq · v² · gap / σ_AN_water   [σ_AN_water = 10.5 mN/m]
```

### Fitting stages planned
- F.1: Transport fit — minimize voltage residual, fit (σ_e, R_contact, V_CE, K_δ) on subset of rows
- F.2: Kinetics fit — minimize FE residual, fit (j0_r, α_c,r) with transport frozen
- F.3: Joint refinement on all rows + 20% holdout

### What needs to be implemented
1. `cell_voltage.jl` — R_series, V_cathode_from_cell functions
2. `hydrodynamics.jl` — delta_lam, We_number, load_bloomquist_data
3. `bloomquist_data.csv` — actual experimental data (user has paper + SI)
4. `run_delta_sweep.jl` — run model at multiple δ values (25, 75, 100, 125 μm) to build interpolation table
5. `plot_fitting.py` — fitting pipeline and comparison plots
6. Possibly `fixed_j_solver.jl` — interpolate (δ, V) → (j, FE) for given j_target

**The user has the Bloomquist paper, SI, and data table and will provide actual data.**

---

## 8. Important Parameters (`params.jl` values)

```julia
F       = 96485.332  C/mol
R_gas   = 8.314463   J/(mol·K)
T       = 298.15     K

K_w     = 1.0e-8     (mol/m³)²   # pKw=14 in M units
K_a2    = 6.3e-5     mol/m³      # pKa2=7.20
K_a3    = 4.5e-10    mol/m³      # pKa3=12.35
k1f     = 1.4        mol/(m³·s)
k1r     = 1.4e8      m³/(mol·s)
k2f     = 1.0e5      m³/(mol·s)
k2r     = 15.87      s⁻¹
k3f     = 2.0e3      m³/(mol·s)
k3r     ≈ 44444      s⁻¹

M_AN    = 0.05306    kg/mol
RHO_AN  = 806.0      kg/m³
m_AN    = 11.59      (partition coeff, Suwanvaipattana 2017)
C_AN_SAT        ≈ 1310   mol/m³
MOLAR_DENSITY_AN ≈ 15191  mol/m³
EPS_ORG_SAT      ≈ 0.0862

C_P_total  = 500.0   mol/m³  (0.5 M Na₃PO₄)
c_Na_input = 1520.0  mol/m³  (3×500 Na⁺ + 20 TBA⁺)
```

---

## 9. How to Run

```
# Stage 1 (ε=0.02, δ=50 μm)
julia an_ehd/run_stage1.jl

# Stage 2 (arithmetic D_mix, ε sweep)
julia an_ehd/run_stage2.jl

# Stage 2m (m_i-corrected D_eff)
julia an_ehd/run_stage2m.jl

# Comparison plots (Python)
cd an_ehd
python plot_stage2_comparison.py
```

All output goes to `an_ehd/output/data/` (CSVs) and `an_ehd/output/plots/` (PNGs).

---

## 10. User Preferences and Constraints

- **Do only what is asked. Do not implement things that weren't requested.**
- Windows WDAC policy blocks many Julia packages from loading (Roots.jl confirmed blocked). Use only Base Julia + LinearAlgebra + SparseArrays + ForwardDiff in solver code.
- Python plotting scripts use: numpy, pandas, matplotlib, glob, re, os — standard only.
- The user will review outputs at each stop gate before proceeding.
- Do not add comments explaining what code does — only add comments for non-obvious WHY (hidden constraints, workarounds).
- Responses should be short and concise.
