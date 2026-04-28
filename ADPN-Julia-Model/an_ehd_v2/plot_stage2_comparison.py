"""
plot_stage2_comparison.py — Arithmetic vs m_i-corrected D_eff comparison.

Loads both Stage 2 (arithmetic D_mix) and Stage 2m (m_i-corrected D_eff) records
and overlays them to visualise the effect of the partition coefficient on:
  1. Peak FE_ADPN vs ε_org
  2. FE_ADPN vs j (Bloomquist comparison)
  3. D_AN,mix vs ε_org
  4. Polarization curves at selected ε_org
  5. AN surface depletion vs j per ε_org
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_DIR  = os.path.join(HERE, "output", "data")
PLOT_DIR  = os.path.join(HERE, "output", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

EPS_SAT = 1310.0 / (806.0 / 0.05306)   # = 0.0862 — Convention A

# Partition coefficients (from params.jl m_partition) for diagnostic curve
m_AN, m_ADPN, m_PN = 11.59, 7.72, 11.59
D_AN_aq, D_AN_org = 2.3e-9, 6.0e-9
D_OH_aq, D_OH_org = 5.27e-9, 0.0

def load_stage(tag):
    """Load all records CSVs for the given stage tag. Returns list of (eps, df)."""
    files = sorted(glob.glob(os.path.join(DATA_DIR, f"{tag}_records_eo*.csv")))
    sweeps = []
    for f in files:
        m = re.search(r"eo(\d+\.\d+)", os.path.basename(f))
        if not m:
            continue
        eps = float(m.group(1))
        df = pd.read_csv(f).sort_values("V", ascending=False).reset_index(drop=True)
        df["j_mA"] = df.j_total * 0.1
        df["j1_mA"] = df.j1 * 0.1
        df["j3_mA"] = df.j3 * 0.1
        sweeps.append((eps, df))
    sweeps.sort(key=lambda x: x[0])
    return sweeps

s2  = load_stage("stage2")
s2m = load_stage("stage2m")
eps_vals = sorted(set([e for e, _ in s2] + [e for e, _ in s2m]))
print(f"Loaded arithmetic sweeps: {[e for e, _ in s2]}")
print(f"Loaded m-corrected sweeps: {[e for e, _ in s2m]}")

# Color map by ε_org
cmap = plt.cm.viridis
colors = {e: cmap(0.15 + 0.75 * i / max(len(eps_vals) - 1, 1))
          for i, e in enumerate(eps_vals)}

# ---------------------------------------------------------------
# 1. Peak FE_ADPN vs ε_org — both formulations on one axis
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5))

def peak_fe(sweeps):
    return np.array([[e, df.FE_ADPN.max(),
                       df.loc[df.FE_ADPN.idxmax(), "V"],
                       df.loc[df.FE_ADPN.idxmax(), "j_mA"]]
                      for e, df in sweeps])

pk_s2  = peak_fe(s2)
pk_s2m = peak_fe(s2m)

ax.plot(pk_s2[:, 0],  pk_s2[:, 1],  "-o", color="C0", ms=8, lw=1.8,
        label="arithmetic D_mix (Stage 2)")
ax.plot(pk_s2m[:, 0], pk_s2m[:, 1], "-s", color="C3", ms=8, lw=1.8,
        label="m_i-corrected D_eff (Stage 2m)")
ax.axvline(EPS_SAT, color="gray", ls="--", lw=0.8, label=f"ε_sat = {EPS_SAT:.4f}")
ax.axhspan(73, 80, alpha=0.15, color="C2", label="Bloomquist target (73–80%)")
ax.set_xlabel(r"$\varepsilon_{\rm org}$")
ax.set_ylabel("Peak FE_ADPN [%]")
ax.set_title("Effect of partition-corrected D_eff on peak FE_ADPN")
ax.set_ylim(0, 100)
ax.legend(fontsize=9, loc="best")
ax.grid(True, ls=":")
# annotate each point with current
for e, fe, V, j in pk_s2m:
    ax.annotate(f" j={j:.0f} mA", (e, fe), fontsize=7, alpha=0.7)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2vs2m_FE_peak.png"), dpi=160)
print("wrote stage2vs2m_FE_peak.png")

# ---------------------------------------------------------------
# 2. FE_ADPN vs j — dual comparison
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
for ax, (sweeps, title) in zip(axes, [(s2, "arithmetic D_mix (Stage 2)"),
                                        (s2m, "m_i-corrected D_eff (Stage 2m)")]):
    for e, df in sweeps:
        regime = "single" if e < EPS_SAT else "two-phase"
        mk = "o" if e < EPS_SAT else "s"
        ax.semilogx(df.j_mA.abs() + 1e-6, df.FE_ADPN, f"-{mk}",
                    color=colors[e], ms=4, lw=1.4,
                    label=f"ε={e:.3f} ({regime})")
    ax.set_xlabel(r"$|j_{\rm total}|$ [mA cm$^{-2}$]")
    ax.set_ylabel("FE_ADPN [%]")
    ax.set_title(title)
    ax.set_ylim(-2, 102)
    ax.axvspan(70, 300, alpha=0.10, color="gray")
    ax.axhspan(73, 80, alpha=0.10, color="C2")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", ls=":")
fig.suptitle("Bloomquist comparison — arithmetic vs m_i-corrected D", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2vs2m_FE_vs_j.png"), dpi=160)
print("wrote stage2vs2m_FE_vs_j.png")

# ---------------------------------------------------------------
# 2b. FE_ADPN vs V — dual comparison
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
for ax, (sweeps, title) in zip(axes, [(s2, "arithmetic D_mix (Stage 2)"),
                                        (s2m, "m_i-corrected D_eff (Stage 2m)")]):
    for e, df in sweeps:
        regime = "single" if e < EPS_SAT else "two-phase"
        mk = "o" if e < EPS_SAT else "s"
        ax.plot(-df.V, df.FE_ADPN, f"-{mk}",
                color=colors[e], ms=4, lw=1.4,
                label=f"ε={e:.3f} ({regime})")
    ax.set_xlabel(r"$-V$ vs SHE  [V]")
    ax.set_ylabel(r"FE$_{\rm ADPN}$ [%]")
    ax.set_title(title)
    ax.set_ylim(-2, 102)
    ax.axhspan(73, 80, alpha=0.10, color="C2", label="Bloomquist 73–80%")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, ls=":")
fig.suptitle(r"FE$_{\rm ADPN}$ vs V — arithmetic vs m_i-corrected D", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2vs2m_FE_vs_V.png"), dpi=160)
print("wrote stage2vs2m_FE_vs_V.png")

# ---------------------------------------------------------------
# 3. D_AN,mix vs ε_org — both formulations
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5))
eps_fine = np.linspace(0, 0.35, 400)
D_AN_arith = np.where(eps_fine < EPS_SAT, D_AN_aq,
                       eps_fine * D_AN_org + (1 - eps_fine) * D_AN_aq)
D_AN_mcorr = np.where(eps_fine < EPS_SAT, D_AN_aq,
                       eps_fine * m_AN * D_AN_org + (1 - eps_fine) * D_AN_aq)
ax.plot(eps_fine, D_AN_arith * 1e9, "C0-", lw=2, label="arithmetic D_AN,mix")
ax.plot(eps_fine, D_AN_mcorr * 1e9, "C3-", lw=2, label="m-corrected D_AN,eff")
ax.axvline(EPS_SAT, color="gray", ls="--", lw=0.8, label=f"ε_sat = {EPS_SAT:.4f}")
ax.set_xlabel(r"$\varepsilon_{\rm org}$")
ax.set_ylabel(r"D$_{\rm AN}$ [10$^{-9}$ m$^2$/s]")
ax.set_title(f"AN diffusivity — arithmetic vs m_i={m_AN}-corrected")
ax.legend(fontsize=9)
ax.grid(True, ls=":")
# annotate factors at sweep points
for e in [0.09, 0.15, 0.25, 0.30]:
    D_a = e * D_AN_org + (1 - e) * D_AN_aq
    D_m = e * m_AN * D_AN_org + (1 - e) * D_AN_aq
    ax.annotate(f" ×{D_m/D_a:.1f}", (e, D_m * 1e9), fontsize=8, color="C3")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2vs2m_D_AN.png"), dpi=160)
print("wrote stage2vs2m_D_AN.png")

# ---------------------------------------------------------------
# 4. Polarization at ε_org = 0.15 — arithmetic vs m-corrected
# ---------------------------------------------------------------
def find_eps(sweeps, target):
    best = min(sweeps, key=lambda x: abs(x[0] - target))
    return best[1] if abs(best[0] - target) < 1e-3 else None

target_eps = 0.15
df_s2  = find_eps(s2, target_eps)
df_s2m = find_eps(s2m, target_eps)

if df_s2 is not None and df_s2m is not None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for df, label, col in [(df_s2,  "arithmetic",    "C0"),
                             (df_s2m, "m-corrected",  "C3")]:
        negV = -df.V
        axes[0].semilogy(negV, df.j1_mA, "-o", color=col, ms=4, label=label)
        axes[1].semilogy(negV, df.j3_mA, "-o", color=col, ms=4, label=label)
        axes[2].plot(    negV, df.FE_ADPN, "-o", color=col, ms=4, label=label)
    axes[0].set(title=f"j_ADPN at ε_org={target_eps}", xlabel="−V [V vs SHE]",
                ylabel=r"j [mA cm$^{-2}$]"); axes[0].legend(); axes[0].grid(True, which="both", ls=":")
    axes[1].set(title=f"j_HER at ε_org={target_eps}", xlabel="−V [V vs SHE]",
                ylabel=r"j [mA cm$^{-2}$]"); axes[1].legend(); axes[1].grid(True, which="both", ls=":")
    axes[2].set(title=f"FE_ADPN at ε_org={target_eps}", xlabel="−V [V vs SHE]",
                ylabel="FE_ADPN [%]", ylim=(-2, 102))
    axes[2].legend(); axes[2].grid(True, ls=":")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f"stage2vs2m_polar_eps{target_eps}.png"), dpi=160)
    print(f"wrote stage2vs2m_polar_eps{target_eps}.png")

# ---------------------------------------------------------------
# 5. AN depletion — does higher D_eff keep surface AN replenished?
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for ax, (sweeps, title) in zip(axes, [(s2, "arithmetic"),
                                        (s2m, "m-corrected")]):
    for e, df in sweeps:
        mk = "o" if e < EPS_SAT else "s"
        ax.semilogx(df.j_mA.abs() + 1e-6, df.AN_depletion, f"-{mk}",
                    color=colors[e], ms=3, lw=1.3, label=f"ε={e:.3f}")
    ax.set_xlabel(r"$|j_{\rm total}|$ [mA cm$^{-2}$]")
    ax.set_ylabel(r"$c_{\rm AN}(0)/c_{\rm AN,bulk}$")
    ax.set_title(title)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", ls=":")
fig.suptitle("AN surface depletion — does m-corrected D delay depletion?", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2vs2m_AN_depletion.png"), dpi=160)
print("wrote stage2vs2m_AN_depletion.png")

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
out_txt = os.path.join(PLOT_DIR, "stage2vs2m_summary.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("Stage 2 vs Stage 2m — arithmetic D_mix vs m_i-corrected D_eff\n")
    f.write("=" * 78 + "\n")
    f.write(f"ε_sat = {EPS_SAT:.4f}\n")
    f.write(f"m_AN = {m_AN}, m_ADPN = {m_ADPN}, m_PN = {m_PN}\n\n")

    f.write(f"{'ε_org':>7} {'regime':>12} "
            f"{'peak FE (arith)':>16} {'peak FE (m-corr)':>17} {'ΔFE':>8}\n")
    f.write("-" * 78 + "\n")
    pk2  = {e: (fe, V, j) for (e, df) in s2
            for fe, V, j in [(df.FE_ADPN.max(),
                              df.loc[df.FE_ADPN.idxmax(), "V"],
                              df.loc[df.FE_ADPN.idxmax(), "j_mA"])]}
    pk2m = {e: (fe, V, j) for (e, df) in s2m
            for fe, V, j in [(df.FE_ADPN.max(),
                              df.loc[df.FE_ADPN.idxmax(), "V"],
                              df.loc[df.FE_ADPN.idxmax(), "j_mA"])]}
    for e in sorted(set(pk2) | set(pk2m)):
        regime = "single" if e < EPS_SAT else "two-phase"
        fe_a = pk2.get(e,  (None,))[0]
        fe_m = pk2m.get(e, (None,))[0]
        if fe_a is None:  fe_a_s = "    —   "
        else:             fe_a_s = f"    {fe_a:>6.2f}% "
        if fe_m is None:  fe_m_s = "    —   "
        else:             fe_m_s = f"    {fe_m:>6.2f}% "
        delta = (fe_m - fe_a) if (fe_a is not None and fe_m is not None) else float("nan")
        f.write(f"{e:>7.3f} {regime:>12} {fe_a_s:>16} {fe_m_s:>17} {delta:>+7.2f}\n")
    f.write(f"\nBloomquist target: FE_ADPN = 73-80% at ε_org ~ 0.15, j = 70-300 mA/cm²\n")
    f.write(f"D_AN boost at ε=0.15: arithmetic = 2.86e-9 → m-corrected = 12.4e-9 (5.4x)\n")
print(f"wrote {out_txt}")

plt.close("all")
print("\nComparison plots complete →", PLOT_DIR)
