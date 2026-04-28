# v7 implementation handoff

**For:** the collaborator's Claude session picking up v7 implementation.
**From:** the previous session (2026-04-27 / 28). The v6.x fit shipped, n_ADN/n_PN are now fit parameters, n_e_TCH was confirmed as 2 from the Bloomquist SI, and the v7 spec is written. Your job is the code.

---

## Read first

Read these in order before touching code:

1. **[Guide Docs/ADPN_EHD_Implementation_Guide_v7.md](../Guide%20Docs/ADPN_EHD_Implementation_Guide_v7.md)** — the v7 spec. Self-contained for what changes from v6.
2. **[Guide Docs/ADPN_EHD_Implementation_Guide_v6.md](../Guide%20Docs/ADPN_EHD_Implementation_Guide_v6.md)** — canonical for sections v7 doesn't touch (governing equations, transport, buffer chemistry, BCs, numerics, cell voltage, hydrodynamics).
3. **[Guide Docs/ADPN_EHD_Implementation_Guide_CHANGELOG_V5toV6.md](../Guide%20Docs/ADPN_EHD_Implementation_Guide_CHANGELOG_V5toV6.md)** — design context. §9 has the v6.x first-fit results; the n_ADN @ LB diagnostic is what motivated adding TCH.
4. **[CONTEXT_TRANSFER.md](CONTEXT_TRANSFER.md)** — older session-handoff notes from before v6.

Skim Bloomquist's paper + SI in this folder for the FE_TCH column definition (see open question §1 below).

---

## Where the code lives

- **`an_ehd_v2/`** — v7 model code. Source files copied from `an_ehd/` as the starting point. Empty `output/` ready.
- **`an_ehd/`** — v6.x baseline. Don't edit. `run_stage4.jl` and `run_stage4v2.jl` are gated with `error()` to prevent overwriting frozen output.
- **`Experimental_data/`** at project root (moved from `an_ehd/Experimental_data/`). Both `an_ehd/` and `an_ehd_v2/` should read from this shared path. **Path update needed in v7:** any file referencing `joinpath(@__DIR__, "Experimental_data", ...)` must become `joinpath(@__DIR__, "..", "Experimental_data", ...)`.

---

## Open questions to resolve before code (or run with TODO defaults)

These are the `<TODO: ...>` flags scattered through the v7 guide. If the user provides values, swap them in; otherwise the defaults are reasonable starting points but should be flagged as approximations in commit messages.

1. **TCH FE column definition in Bloomquist SI.** Verify what `n_e_TCH` and `MW_TCH` Bloomquist used to compute the `FE_TCH_pct` column in `Experimental_data/bloomquist_data.csv`. The user reports the SI says `n_e_TCH = 2`. If Bloomquist computed FE_TCH with MW = 175.23 (saturated C₉H₁₃N₃), but `n_e = 2` mass-balance gives MW = 161.20 (C₉H₁₁N₃), there is a ~1.087× systematic correction needed on FE_TCH residuals. **Highest priority cross-check.**
2. **TCH partition coefficient `m_TCH`** (`an_ehd_v2/params.jl`). v7 default is 1.5 (estimate); this is the single most uncertain transport input and the user said they'd ask the lab.
3. **TCH diffusivities** `D_TCH,aq`, `D_TCH,org`. v7 defaults: 7.0×10⁻¹⁰, 1.2×10⁻⁹ m²/s. Wilke-Chang scaling from ADPN — fine starting point.
4. **TCH standard potential `E°_TCH`**. v7 default: −1.30 V vs SHE (same as ADPN/PN). Bui Lab may have a TCH-specific value.

---

## Implementation tasks (ordered)

Estimated effort and rough sequencing. Each task should be its own commit.

### Phase 1 — Core 9-species model (~3 hours)

1. **`an_ehd_v2/params.jl`** — add TCH constants, extend `D_aq`, `D_org`, `m_partition`, `z_species`, `c_bulk` to length 9. Bump `n_species = 9`, `JAC_BLOCK = 10`, `JAC_HALFBW = 19`, `n_colors = 39`. Add `j0_TCH`, `alpha_c_TCH`, `n_TCH_default = 3.0`, `E0_TCH = −1.30`, `nE_TCH = 2`. Update HER defaults to v6.x converged values (`j0_3 = 2.666e-5`, `alpha_c3 = 0.390`).
2. **`an_ehd_v2/kinetics.jl`** — extend `KIN_OVERRIDE` Ref to `j0::NTuple{4}, ac::NTuple{4}, n::NTuple{3}`. `tafel_currents` returns 4-tuple `(j_1, j_2, j_3, j_4)`. Add `j_4 = j0_TCH · (cA/c_ref)^n_TCH · exp(−α_TCH · F · η_TCH / RT)`.
3. **`an_ehd_v2/diffusivity.jl`** — verify `D_mix` arrays handle length 9 (probably already parametric over `Params.n_species`). Audit for hardcoded 8.
4. **`an_ehd_v2/chemistry.jl`** — `make_initial_guess` and `bulk_concentration` populate index 9 (`c_TCH = 0`).
5. **`an_ehd_v2/transport.jl`** — audit the Scharfetter-Gummel loop for `1:n_species` parametrisation. Probably no change needed, but verify.
6. **`an_ehd_v2/assembly.jl`** — DOF indexing `9·(ix−1)+k → 10·(ix−1)+k`. Update Faradaic flux BCs per v7 §5.1:
   - `J_OH(0) = +(j_1 + j_2 + j_3 + j_4) / F`
   - `J_AN(0) = −(2 j_1 + j_2 + 3 j_4) / (2 F)`
   - `J_TCH(0) = +j_4 / (2 F)` (NEW)
   - `J_ADPN`, `J_PN` unchanged.
   Audit for hardcoded `9` or `8` anywhere.
7. **Smoke test.** Run `an_ehd_v2/test_smoke.jl` to verify Newton converges at V = −1.0 V from cold IC. Fix any issues here before moving on. **Do not skip this gate** — it's the cheapest signal that the 9-species DOF layout is internally consistent.

### Phase 2 — Stage 1 / 2 / 2m / 3 caches (~30 min wake compute)

8. **`an_ehd_v2/run_stage1.jl`** — re-run; verify convergence and that the v7 baseline FE_ADN / FE_PN / FE_TCH match v6 baseline within ~1 pp (v7 model with default kinetics should reproduce v5/v6 behaviour where AN orders are 2/1/3 and TCH adds a small new sink).
9. **`an_ehd_v2/run_stage2.jl`, `run_stage2m.jl`** — re-run ε_org sweeps. Plots in `an_ehd_v2/output/stage2*/plots/` should look like v6's plots with one extra line (TCH branch).
10. **`an_ehd_v2/run_stage3.jl`** — populate the warm-start cache for v7. Outputs to `an_ehd_v2/output/cache/`. Required before Stage 4v3.

### Phase 3 — Stage 4v3 fit (~75 min wake compute)

11. **`an_ehd_v2/fixed_j_solver.jl`** — `solve_at_j` accepts `j0::NTuple{4}`, `alpha_c::NTuple{4}`, `n_orders::NTuple{3}`. Default `n_orders = (2.0, 1.0, 3.0)`.
12. **`an_ehd_v2/fit_kinetics.jl`** — `N_THETA = 9`. New layout per v7 §8.1. HER frozen via `J0_3_FROZEN = 2.666e-5`, `ALPHA_C3_FROZEN = 0.390` constants near the top of the module; injected when `theta_to_physical` builds the 4-tuples passed to `solve_at_j`. Residuals: `F` length is `3 · length(sel)` with FE_ADN, FE_PN, FE_TCH stacked per row.
13. **`an_ehd_v2/run_stage4v3.jl`** — new file mirroring `run_stage4v2.jl`. Outputs to `an_ehd_v2/output/stage4v3/`. Saves all 9 fit params plus the frozen HER values to `stage4a_fitted_theta.txt`. Update path to bloomquist_data.csv per the moved-data note above.
14. **Run.** `caffeinate -dimsu -w <julia_pid>` and stay on AC. Surface the LM iterations as they land. Final scoreboard reports gate pass/fail for FE_ADN < 8 pp, FE_PN < 5 pp, FE_TCH < 4 pp, Extended ADN < 12 pp, Holdout ADN < 15 pp.

### Phase 4 — Diagnostic plots (~30 min)

15. **`an_ehd_v2/plot_stage4v3_parity.py`** — 4-panel composite (ADN, PN, TCH, V_cell). Update `DATA_FN` to `an_ehd_v2/output/stage4v3/data/stage4_diagnostic.csv`.
16. **`an_ehd_v2/plot_stage4v3_residuals.py`** — adds FE_TCH residual vs (j, ε_org), faceted by gap.
17. **`an_ehd_v2/plot_stage4v3_3d_surfaces.py`** — adds 4th panel for FE_TCH regime.
18. **`an_ehd_v2/analyze_stage4.jl`** — extend to read 9-key fitted-theta and emit `stage4v3_diagnostic.csv` with FE_TCH column.

### Phase 5 — Documentation (~15 min)

19. Append §10 ("Third fit — TCH-aware, HER-frozen, 2026-04-?? results") to `Guide Docs/ADPN_EHD_Implementation_Guide_CHANGELOG_V5toV6.md` (or split into a new `CHANGELOG_V6toV7.md` if it gets long). Mirror the structure of §9: code changes, run summary, fitted parameters, decision-gate scoreboard, loss trajectory, three findings, output artefacts, next-step priority.

---

## Things to flag if you hit them

- **n_ADN still pinning at LB = 1.0 after TCH lands.** That would point at Langmuir-Hinshelwood coverage saturation rather than the missing-TCH hypothesis. Don't relax the bound mid-fit; finish v3, document the result, then design v4 with `n_ADN ∈ [0.5, 3.0]` as a follow-up diagnostic.
- **n_TCH pinning at UB = 3.0.** That would mean the trimer-molecularity reaction order is the fit's preferred answer — physically reasonable, but check the residuals look random. If pinned, treat similarly to v6.x's bound-pinning diagnostic and consider relaxing the upper bound on a v4 fit.
- **Negative or non-finite FE_TCH values** in the diagnostic CSV. Probably a sign of the c_AN^n_TCH branch hitting numerical noise when c_AN(surface) is very small. The `max(c_AN_surface, 0)` guard in `tafel_currents` should protect, but verify.
- **HER residuals showing systematic bias vs j or ε_org.** If the frozen v6.x HER values fight a meaningful fraction of rows, thaw them back into the fit (return to N_THETA = 11) and re-run as v3b. Document in the changelog.
- **MW_TCH discrepancy.** If you confirm Bloomquist used MW = 175.23 (saturated C₉H₁₃N₃) in their FE_TCH column despite the SI saying `n_e = 2`, multiply model FE_TCH by 175.23/161.20 = 1.087 *before* computing residuals. Flag this in the changelog as a caveat; the residual is otherwise biased.

---

## Open communication channels

- The user (Sebobbit) reviews each stop gate. After Phase 1 smoke test, after Phase 2 Stage 1 baseline, and after Phase 3 Stage 4v3 results — pause and ask for review.
- Caffeinate gotcha: the `-s` flag is documented as AC-only on macOS. Lid closure on battery suspends the Julia process even with `caffeinate -dimsu -w <PID>` running. Tell the user explicitly to stay on AC during the Stage 4v3 fit.

---

## Status as of handoff (2026-04-28)

- v6 fit baseline frozen in `an_ehd/output/stage4/` (CHANGELOG §8).
- v6.x fit baseline frozen in `an_ehd/output/stage4v2/` (CHANGELOG §9). Loss 7.51e+03; Core RMSE FE_ADN 11.32 pp / FE_PN 5.33 pp; Holdout 29.33 pp; n_ADN @ LB.
- v6.x patches merged into `an_ehd/`: kinetics override carries `n=(n_ADN, n_PN)` 2-tuple; `tafel_currents` uses `(c/c_ref)^n_r` form; fit_kinetics N_THETA = 8.
- `n_1` / `n_2` renamed to `n_ADN` / `n_PN` across kinetics, fit_kinetics, run_stage4v2.jl.
- `an_ehd/run_stage4.jl` gated with `error()` to protect v6 baseline.
- `an_ehd_v2/` skeleton populated with copies of `an_ehd/*.jl` and `an_ehd/*.py`. **Output empty.** Tested: nothing yet.
- Bloomquist SI confirms `n_e_TCH = 2` per the user. MW / molecular formula still needs cross-check against the FE_TCH column derivation (see open question §1 above).

Good luck. The v6.x → v7 step is mechanically straightforward but has many small places to introduce bandwidth or DOF-index bugs, so favour small commits and run `test_smoke.jl` after each module change.
