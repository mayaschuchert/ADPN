# ADPN EHD Model — Running Changelog

**Bui Lab | NYU Tandon School of Engineering**

Tracks changes from `ADPN_EHD_Implementation_Guide_v3.md` (work-in-progress) to the current working version `v4.md`.

---

## v3 → v4 Changes

### Scope and Temperature

| Item | v3 | v4 |
|------|----|----|
| Operating temperature | Not stated | Pinned to T = 298.15 K (25 °C), added to scope line, parameter table §9.1, sg_flux signature, and all K value comments |
| T as Julia constant | Absent | `const T = 298.15` added to params.jl block |

All equilibrium constants (K_w, K_a2, K_a3), diffusivities, and kinetic expressions are now explicitly at 25 °C.

---

### Species and Degrees of Freedom (§2)

| Item | v3 | v4 |
|------|----|----|
| Species count | 10 (including TBA⁺) | 9 (TBA⁺ removed) |
| TBA⁺ treatment | Listed as independent species (z = +1, fixed c = 20 mol m⁻³) | Lumped into Na⁺ pool; no separate tracking |
| Effective Na⁺ input | 3 × 500 = 1500 mol m⁻³ | 3 × 500 + 20 = **1520 mol m⁻³** (absorbs TBA⁺ contribution) |
| DOFs per cell | 10 (9 log-concentrations + φ_l) | **9** (8 log-concentrations + φ_l) |
| Total unknowns (N = 100) | — | **900** |
| Electroneutrality formula | Included TBA⁺ term | No c_TBA term; c_Na = c_OH + c_H₂PO₄ + 2·c_HPO₄ + 3·c_PO₄ − c_H |

**Rationale for lumping TBA⁺ into Na⁺:** Both carry z = +1. Na⁺ is not transported — it is recovered from electroneutrality at each cell — so no mobility mismatch enters the model. Lumping removes a DOF without any physical approximation beyond the Kohlrausch assumption already made by the electroneutrality closure.

---

### D_aq / D_org Arrays (§4)

| Item | v3 | v4 |
|------|----|----|
| Array length | 9 (included TBA⁺ entry) | **8** (H⁺, OH⁻, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻, AN, ADPN, PN) |
| TBA⁺ diffusivity | Placeholder entry present | Removed |

---

### Buffer Chemistry (§6)

| Item | v3 | v4 |
|------|----|----|
| buffer_sources! index range | R[6:9] (included TBA⁺ slot) | **R[6:8]** (AN, ADPN, PN only) |
| Bulk equilibrium solver | Absent | **New §6.4** — `solve_phosphate_equilibrium()` via bisection on charge balance |
| Initial guess builder | Absent | **New** — `make_initial_guess()`: flat profile at equilibrium concentrations, φ_l = 0 |
| Verification assertions | Absent | Added: K_w/K_a2/K_a3 product checks; zero buffer residual check |
| Expected equilibrium values | Absent | Added: table (pH ≈ 13.0, c_PO₄ ≈ 410–430 mol m⁻³, etc.) |

---

### Boundary Conditions (§7)

| Item | v3 | v4 |
|------|----|----|
| Bulk BC source | Not specified | Explicitly references `solve_phosphate_equilibrium()` |
| φ_l(δ) = 0 role | Not explained | Added note: serves as **gauge fix** for current conservation null direction |
| φ_l(0) | Not discussed | Clarified: not prescribed — determined by current conservation equation |

---

### Numerical Methods (§10)

| Item | v3 | v4 |
|------|----|----|
| N_mesh | 200 | **100** |
| Mesh definition | Fixed absolute dx_min, dx_max | **Stretch-factor parameterisation**: s = dx_max/dx_min = 10; dx_min computed from geometric series sum as function of δ and N |
| Mesh function | None | **`make_mesh(N, delta; stretch=10.0)`** — fills [0, δ] exactly for any δ |
| Jacobian strategy | Dense FD (900+ columns) | **Banded FD with column coloring**: 19 perturbation vectors for bandwidth-18 block-tridiagonal structure |
| Jacobian solver | Dense LU | **Sparse/banded LU** (SparseArrays.jl + KLU.jl) |
| Jacobian cost | ~900 residual evaluations / step | ~19 residual evaluations / step (~47× speedup) |
| φ_l step clamp | 0.1 V | **0.015 V** |
| Continuation strategy | Pseudo-arclength continuation (PAC) | **Simple Newton continuation** (PAC removed as primary strategy) |
| Continuation step control | PAC adaptive arclength | ds = 0.05 V; increase ×1.5 if ≤ 5 Newton iterations; halve and retry if fails; ds_min = 0.005 V |
| Residual form note | Standalone callout | Integrated into §10.5 with updated dx-ratio note (s = 10, not ~100) |

**Rationale for N = 100 (was 200):** The Scharfetter–Gummel scheme is exact in 1D for uniform fields between nodes; resolution need only resolve concentration boundary layers near x = 0. At δ = 50 μm with stretch = 10, dx_min ≈ 0.1 μm, which resolves the ~1 μm diffusion sub-layer near the electrode. Halving N cuts Jacobian assembly cost by 4× (banded: O(N)) and LU cost by ~2× (banded LU: O(N × b²)).

**Rationale for banded FD (was dense):** The block-tridiagonal Jacobian has half-bandwidth b = 9. Column-grouping with 2b+1 = 19 colors reduces residual evaluations per Jacobian from O(N) to a fixed 19, independent of mesh size. Combined with sparse LU, this brings per-step cost from ~60 ms to ~0.3 ms at N = 100.

**Rationale for φ_l clamp 0.015 V (was 0.1 V):** Physical φ_l variation across the diffusion layer is a few mV at most. A 0.1 V clamp permits potential steps far exceeding the physical range, which couple into SG fluxes via exp(z·F·Δφ/RT) and destabilise concentration DOFs.

**Rationale for Newton continuation (was PAC):** PAC is warranted when the solution curve has fold bifurcations (turning points in the continuation parameter). The 1D NP diffusion layer is a monotone system — current increases monotonically with overpotential, with no bistability. Simple Newton continuation with adaptive step size is sufficient and avoids the 2× system-size augmentation that PAC requires.

---

### Bootstrap Protocol (§12)

| Item | v3 | v4 |
|------|----|----|
| α_buf ramp | Vague ("increase from 0 → 1") | **10 uniform steps**, Δα = 0.1, at V = −1.0 V |
| α_kin ramp | "geometrically" (vague) | **Geometric ×2 per step**: 10⁻⁶, 2×10⁻⁶, …, 1.0 (~20 steps) |
| V sweep method | PAC sweep | Newton continuation sweep (§10.4) |

---

### Module Structure (§14)

| Item | v3 | v4 |
|------|----|----|
| mesh.jl description | "1D graded mesh" | `make_mesh(N, delta; stretch)` noted explicitly |
| solver.jl description | `newton_solve!, pac_solve` | `newton_solve!` (banded FD, sparse LU), `newton_continuation`; PAC removed |
| chemistry.jl description | `buffer_sources!` only | Added `solve_phosphate_equilibrium`, `make_initial_guess` |
| Sparse solver note | Absent | Added: "no dense Jacobian"; KLU.jl noted |

---

### Pitfalls (§15)

New rows added in v4:

| New Pitfall | Symptom | Fix |
|-------------|---------|-----|
| Cold start (no equilibrium IC) | Newton diverges immediately | Use `solve_phosphate_equilibrium` + `make_initial_guess` |
| Divided FV form on graded mesh | Jacobian ill-conditioned | Use integrated form J_L − J_R + S·dx |
| Missing φ_l gauge fix | Jacobian rank-deficient | Set φ_l(δ) = 0 (Dirichlet) |
| j₀ in mA/cm² passed directly | Rates off by ×10 | Multiply by 10: A/m² = mA/cm² × 10 |
| Dense FD Jacobian | Full sweep >10 min | Banded FD, 19-column grouping, sparse LU |
| Fixed absolute dx_min/dx_max | Domain length wrong across δ sweep | Stretch-factor mesh via make_mesh |
| φ_l clamp at 0.1 V | Over-large steps destabilise c_i | Clamp to 0.015 V |
| Implementing PAC prematurely | Unnecessary complexity | Use simple Newton continuation first |

---

### Other Additions (§3, §9)

| Item | v3 | v4 |
|------|----|----|
| §3.2 current conservation | Equation present, no explanation | Added note: null direction in φ_l, gauge fix at φ_l(δ) = 0 required for non-singular Jacobian |
| §9.1 Physical Constants | Absent | Added: F, R, T table |
| §9.3 electrolyte | TBA-OH listed without note | Added "lumped into Na⁺" annotation |
| §9.3 EDTA | No note | Added "not tracked (negligible vs. phosphate)" |

---

## v4 Concrete Implementation Pass

After the v3 → v4 migration, the guide still listed component pieces (sg_flux, buffer_sources!, make_initial_guess) but never showed how they assemble into a runnable model. This pass added the executable spine.

### Continuation Parameters Defined (§10.6, new)

α_buf and α_kin had been mentioned in §12 but never defined. Added explicit definitions:

| Parameter | Scales | At α = 0 | At α = 1 |
|-----------|--------|----------|----------|
| α_buf | Buffer source terms R_buf,i | Buffer chemistry off (equilibrium IC is exact) | Full buffer kinetics |
| α_kin | Tafel current densities j₁, j₂, j₃ | No Faradaic current; flat profiles preserved | Full kinetics |

`buffer_sources!` and `tafel_currents` updated to take `alpha_buf` / `alpha_kin` as multiplicative scalars (defaults 1.0). `buffer_sources!` in §6.2 also updated to match.

### Full Residual Assembly (§10.7, new)

Added the central `full_residual!` function showing how all pieces compose. The function handles five distinct cases that were previously left to the reader:

| Case | Treatment |
|------|-----------|
| Interior cell, species k | Integrated FV form: F = J_left − J_right + S × dx |
| Interior cell, φ_l | Σ z_i (J_left − J_right) over charged species (k = 1..5) |
| First cell, left face (x = 0) | Faradaic boundary fluxes from `tafel_currents()`, evaluated using c_AN[1] and the *solved* φ_l[1] (implicit coupling) |
| Last cell (x = δ), all DOFs | Dirichlet override: F[k] = u[k] − log(c_bulk_k); F[φ] = φ_l[N] |
| Sign convention at electrode | N_i = +x flux; consumption (AN) → N < 0; production (OH⁻, ADPN, PN) → N > 0 |

Supporting additions:
- **DOF indexing helpers** `conc_dof(ix, k) = 9*(ix−1) + k` and `phi_dof(ix) = 9*ix` to prevent off-by-one errors throughout assembly.
- **`bulk_concentration(k, c_eq, eps_org)`** dispatcher that returns Dirichlet values for each species (equilibrium phosphates from `solve_phosphate_equilibrium`, AN from `c_AN_bulk(eps_org)`, ADPN/PN seed at 10⁻³).
- **`c_AN_bulk(eps_org; eps_sat=0.09)`** function implementing the two regimes from §7.2: dilution scaling below saturation, fixed at C_AN_SAT = 1310 mol m⁻³ above.

### Newton Solver Code (§10.8, new)

Replaced the prose description with concrete `newton_solve!` and `newton_continuation` pseudocode covering:

- Sparse Jacobian construction once per call via `build_banded_pattern(n, 9)` — block-tridiagonal, block size 9, half-bandwidth 9
- Banded FD Jacobian via 19-color grouping (calls `banded_fd_jacobian!` from §10.2)
- Sparse LU solve `du = -(J \ F0)` (defers to KLU.jl or LinearAlgebra.lu on the SparseMatrixCSC)
- Per-DOF-type step clamping inside the iteration loop (5.0 for log-conc, 0.015 V for φ_l)
- L∞ tolerance 10⁻⁴; max 25 iterations
- Continuation step adapter: ds × 1.5 if ≤ 5 iterations, ds / 2 on failure, abort below ds_min = 0.005 V

### Sparsity Pattern Note

Added an explicit note that `build_banded_pattern` builds the SparseMatrixCSC pattern once outside the iteration loop, and `banded_fd_jacobian!` writes into the existing non-zeros — avoiding repeated allocation across Newton iterations.

---

## Next Steps — Process for Continuous Refinement

The guide is complete enough to start implementation. Further paper-level polish has diminishing returns; the remaining latent bugs (sign flips, indexing errors, missing terms) will only surface when running the model.

### Recommended cadence

- **Build-driven refinements (primary).** After each Newton-convergence milestone — buffer ramp, kinetics ramp, first continuation step, first full V sweep — do a short guide pass to record what was actually wrong, and bump this changelog in the same pass.
- **Section-by-section pass** only when a specific topic is uncertain (e.g. "is the current-conservation sign right?"). Use before committing time to build that piece.
- **Whole-document pass** after 5–10 small edits to catch consistency drift (species count, DOF layout, cross-references, subscript conventions).

### First build target

**Stage 1, bootstrap step 1 only:** α_buf = 0, α_kin = 0 at V = −1.0 V vs SHE, ε_org = 0, δ = 50 μm.

With the equilibrium initial guess from `make_initial_guess()` (§6.4) and all source/flux terms disabled, Newton should converge in **zero iterations** — the residual is identically zero by construction. This single milestone tests the full infrastructure (residual assembly, DOF layout, mesh generation, Jacobian sparsity pattern, sparse LU) without any physics. If it fails, the bug is in the infrastructure, not the model.

Only after this passes should the α_buf ramp begin.

### On trial and error

Appropriate for **numerical tuning** (step sizes, ramp schedules, tolerances) — these cannot be predicted up front.

Inappropriate for **structural failures** (divergence, non-finite residuals, rank-deficient Jacobian). When something fails structurally, diagnose *why* before tweaking:

- Print ‖F‖_∞ per iteration and identify which DOF row dominates
- Check electroneutrality Σ z_i c_i at each cell
- Plot the concentration and φ_l profiles at the failed iterate
- Verify Σ z_i N_i at each face is near zero at the current iterate

Random parameter changes on a stiff Newton solver will waste days.

### Tradeoff

Starting to build now means some sections of the current guide will turn out wrong and any polish time spent on them is wasted. The reverse risk — spending another week on a perfect-on-paper guide that still has latent bugs — is worse.

---

*Last updated: April 2026*
