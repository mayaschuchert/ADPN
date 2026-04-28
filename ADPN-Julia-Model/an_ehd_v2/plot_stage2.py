"""
plot_stage2.py — Stage 2 ε_org-sweep diagnostic plots.

Inputs:
  output/data/stage2_records_eo*.csv    (one per ε_org)
  output/data/stage2_profile_eo*.csv    (profiles at selected V)
  output/data/stage2_meta_eo*.txt       (per-sweep metadata incl. cAN_bulk)

Outputs (all PNG in output/plots/stage2_*):
  stage2_polarization_overlay.png       j_r, j_total vs V — curves colored by ε_org
  stage2_FE_vs_V.png                    FE_ADPN, FE_PN, FE_HER vs V per ε_org
  stage2_FE_vs_j.png                    FE_ADPN vs j — key Bloomquist comparison (§19 panel c)
  stage2_FE_peak_vs_eps.png             Peak FE_ADPN vs ε_org (§19 panel d)
  stage2_D_vs_eps.png                   D_AN,mix, D_OH,mix vs ε_org — regime transition (§19 g)
  stage2_AN_depletion_vs_j.png          c_AN(0)/c_bulk vs j per ε_org (§19 panel e)
  stage2_phi_vs_j.png                   φ_ℓ(0) vs j per ε_org (§19 panel f)
  stage2_summary.txt                    Text summary of peaks and regime diagnostics
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

# ---------------------------------------------------------------
# Load Stage 2 records (one CSV per ε_org)
# ---------------------------------------------------------------
rec_files = sorted(glob.glob(os.path.join(DATA_DIR, "stage2_records_eo*.csv")))
if not rec_files:
    raise SystemExit("No Stage 2 records found in " + DATA_DIR)

sweeps = []   # list of (eps_org, DataFrame) sorted by eps_org
for f in rec_files:
    m = re.search(r"eo(\d+\.\d+)", os.path.basename(f))
    if not m:
        continue
    eps = float(m.group(1))
    df = pd.read_csv(f).sort_values("V", ascending=False).reset_index(drop=True)
    df["j_total_mA"] = df.j_total * 0.1
    df["j1_mA"]      = df.j1      * 0.1
    df["j2_mA"]      = df.j2      * 0.1
    df["j3_mA"]      = df.j3      * 0.1
    sweeps.append((eps, df))
sweeps.sort(key=lambda x: x[0])

eps_values = [s[0] for s in sweeps]
print(f"Loaded Stage 2 sweeps at ε_org = {eps_values}")

# Convention A saturation threshold (from guide v5 §4.1)
EPS_SAT = 1310.0 / (806.0 / 0.05306)   # ≈ 0.0862

# ε_org colour ramp — light → dark with ε_org
cmap = plt.cm.viridis
colors = {eps: cmap(0.12 + 0.80 * i / max(len(eps_values) - 1, 1))
          for i, eps in enumerate(eps_values)}

def regime_marker(eps):
    return "o" if eps < EPS_SAT else "s"

# ---------------------------------------------------------------
# 1. Polarization overlay — j_r & j_total vs V for each ε_org
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
for eps, df in sweeps:
    neg_V = -df.V
    mk = regime_marker(eps)
    c = colors[eps]
    axes[0].semilogy(neg_V, df.j1_mA, f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
    axes[1].semilogy(neg_V, df.j3_mA, f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
    axes[2].semilogy(neg_V, df.j_total_mA, f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
for i, (title, ylabel) in enumerate([("(a) j_ADPN vs V", r"$j_{\rm ADPN}$ [mA cm$^{-2}$]"),
                                      ("(b) j_HER vs V",  r"$j_{\rm HER}$ [mA cm$^{-2}$]"),
                                      ("(c) j_total vs V", r"$j_{\rm total}$ [mA cm$^{-2}$]")]):
    axes[i].set_xlabel(r"$-V$ vs SHE  [V]")
    axes[i].set_ylabel(ylabel)
    axes[i].set_title(title)
    axes[i].legend(loc="lower right", fontsize=8, ncol=2)
    axes[i].grid(True, which="both", ls=":")
fig.suptitle(r"Stage 2 — polarization curves by $\varepsilon_{\rm org}$  (circle = single-phase, square = two-phase)",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_polarization_overlay.png"), dpi=160)
print("wrote stage2_polarization_overlay.png")

# ---------------------------------------------------------------
# 2. FE vs V — 3 panels (one per reaction), curves by ε_org
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
for eps, df in sweeps:
    neg_V = -df.V
    mk = regime_marker(eps)
    c = colors[eps]
    axes[0].plot(neg_V, df.FE_ADPN, f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
    axes[1].plot(neg_V, df.FE_PN,   f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
    axes[2].plot(neg_V, df.FE_HER,  f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
for i, title in enumerate(["FE_ADPN vs V", "FE_PN vs V", "FE_HER vs V"]):
    axes[i].set_xlabel(r"$-V$ vs SHE  [V]")
    axes[i].set_ylabel("FE [%]")
    axes[i].set_ylim(-2, 102)
    axes[i].set_title(title)
    axes[i].legend(loc="best", fontsize=8, ncol=2)
    axes[i].grid(True, ls=":")
fig.suptitle("Stage 2 — Faradaic efficiencies by ε_org", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_FE_vs_V.png"), dpi=160)
print("wrote stage2_FE_vs_V.png")

# ---------------------------------------------------------------
# 3. FE_ADPN vs j — key Bloomquist comparison (§19 panel c)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5.5))
for eps, df in sweeps:
    mk = regime_marker(eps)
    c = colors[eps]
    regime = "single" if eps < EPS_SAT else "two-phase"
    ax.semilogx(df.j_total_mA.abs(), df.FE_ADPN,
                f"-{mk}", color=c, ms=4, lw=1.4,
                label=f"ε={eps:.3f} ({regime})")
ax.set_xlabel(r"$|j_{\rm total}|$  [mA cm$^{-2}$]")
ax.set_ylabel(r"FE$_{\rm ADPN}$  [%]")
ax.set_ylim(-2, 102)
ax.set_title(r"FE$_{\rm ADPN}$ vs current density — key Bloomquist comparison")
ax.axvspan(70, 300, alpha=0.10, color="gray", label="Bloomquist j-range")
ax.legend(loc="best", fontsize=9)
ax.grid(True, which="both", ls=":")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_FE_vs_j.png"), dpi=160)
print("wrote stage2_FE_vs_j.png")

# ---------------------------------------------------------------
# 4. Peak FE_ADPN vs ε_org (§19 panel d)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.8))
peak_FE_ADPN = [(eps, df.FE_ADPN.max(),
                 df.loc[df.FE_ADPN.idxmax(), "V"],
                 df.loc[df.FE_ADPN.idxmax(), "j_total_mA"])
                for eps, df in sweeps]
eps_arr   = np.array([p[0] for p in peak_FE_ADPN])
peak_arr  = np.array([p[1] for p in peak_FE_ADPN])
V_at_peak = np.array([p[2] for p in peak_FE_ADPN])
j_at_peak = np.array([p[3] for p in peak_FE_ADPN])

ax.plot(eps_arr, peak_arr, "-o", color="C0", ms=7, lw=1.6,
        label="peak FE_ADPN (model)")
ax.axvline(EPS_SAT, color="gray", ls="--", lw=0.8,
           label=f"ε_sat = {EPS_SAT:.4f}")
ax.set_xlabel(r"$\varepsilon_{\rm org}$")
ax.set_ylabel("Peak FE_ADPN [%]")
ax.set_ylim(0, max(80, peak_arr.max() * 1.15))
ax.set_title("Peak FE_ADPN vs organic loading")
ax.legend(fontsize=9)
ax.grid(True, ls=":")
# annotate with V and j at each peak
for eps, fe_peak, V_p, j_p in peak_FE_ADPN:
    ax.annotate(f" V={V_p:.2f}\n j={j_p:.1f} mA",
                (eps, fe_peak), fontsize=7, alpha=0.7)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_FE_peak_vs_eps.png"), dpi=160)
print("wrote stage2_FE_peak_vs_eps.png")

# ---------------------------------------------------------------
# 5. D_AN,mix and D_OH,mix vs ε_org — regime transition (§19 g)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))

# Compute D_mix theoretically for a fine grid + sweep points
D_AN_aq, D_AN_org = 2.3e-9, 6.0e-9
D_OH_aq, D_OH_org = 5.27e-9, 0.0
eps_fine = np.linspace(0, 0.35, 400)
D_AN_line = np.where(eps_fine < EPS_SAT, D_AN_aq,
                     eps_fine * D_AN_org + (1 - eps_fine) * D_AN_aq)
D_OH_line = np.where(eps_fine < EPS_SAT, D_OH_aq,
                     eps_fine * D_OH_org + (1 - eps_fine) * D_OH_aq)
ax.plot(eps_fine, D_AN_line * 1e9, "C0-", lw=1.8, label="D_AN,mix (model)")
ax.plot(eps_fine, D_OH_line * 1e9, "C1-", lw=1.8, label="D_OH,mix (model)")
ax.plot(eps_arr, [D_AN_aq * 1e9 if e < EPS_SAT
                  else (e * D_AN_org + (1 - e) * D_AN_aq) * 1e9
                  for e in eps_arr], "C0o", ms=9, label="sweep points")
ax.plot(eps_arr, [D_OH_aq * 1e9 if e < EPS_SAT
                  else (e * D_OH_org + (1 - e) * D_OH_aq) * 1e9
                  for e in eps_arr], "C1o", ms=9)
ax.axvline(EPS_SAT, color="gray", ls="--", lw=0.8,
           label=f"ε_sat = {EPS_SAT:.4f}")
ax.set_xlabel(r"$\varepsilon_{\rm org}$")
ax.set_ylabel(r"D$_{\rm mix}$  [10$^{-9}$ m$^2$/s]")
ax.set_title("Regime transition — D_AN,mix (↑) and D_OH,mix (↓) at ε_sat")
ax.legend(fontsize=9)
ax.grid(True, ls=":")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_D_vs_eps.png"), dpi=160)
print("wrote stage2_D_vs_eps.png")

# ---------------------------------------------------------------
# 6. AN depletion vs j — (§19 panel e)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
for eps, df in sweeps:
    mk = regime_marker(eps)
    c = colors[eps]
    ax.semilogx(df.j_total_mA.abs() + 1e-6, df.AN_depletion,
                f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
ax.set_xlabel(r"$|j_{\rm total}|$ [mA cm$^{-2}$]")
ax.set_ylabel(r"$c_{\rm AN}(0)\,/\,c_{\rm AN, bulk}$")
ax.set_ylim(-0.05, 1.05)
ax.set_title("AN depletion at electrode vs current")
ax.legend(loc="best", fontsize=8, ncol=2)
ax.grid(True, which="both", ls=":")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_AN_depletion_vs_j.png"), dpi=160)
print("wrote stage2_AN_depletion_vs_j.png")

# ---------------------------------------------------------------
# 7. φ_ℓ(0) vs j — ohmic penalty (§19 panel f)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
for eps, df in sweeps:
    mk = regime_marker(eps)
    c = colors[eps]
    ax.semilogx(df.j_total_mA.abs() + 1e-6, df.phi_l_surface * 1000,
                f"-{mk}", color=c, ms=3, label=f"ε={eps:.3f}")
ax.set_xlabel(r"$|j_{\rm total}|$ [mA cm$^{-2}$]")
ax.set_ylabel(r"$\phi_\ell(0)$ [mV]")
ax.set_title("Electrolyte potential at surface — ohmic penalty")
ax.legend(loc="best", fontsize=8, ncol=2)
ax.grid(True, which="both", ls=":")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage2_phi_vs_j.png"), dpi=160)
print("wrote stage2_phi_vs_j.png")

# ---------------------------------------------------------------
# 8. Summary text
# ---------------------------------------------------------------
out_txt = os.path.join(PLOT_DIR, "stage2_summary.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("Stage 2 summary — ε_org sweep at δ = 50 μm\n")
    f.write("=" * 78 + "\n")
    f.write(f"ε_sat = {EPS_SAT:.4f}\n")
    f.write(f"Number of sweeps: {len(sweeps)}\n\n")
    f.write(f"{'ε_org':>7} {'regime':>12} {'N_pts':>6} {'peak FE_ADPN':>14} "
            f"{'@ V':>8} {'@ j [mA/cm²]':>16}\n")
    f.write("-" * 78 + "\n")
    for (eps, df), (_, fe_peak, V_p, j_p) in zip(sweeps, peak_FE_ADPN):
        regime = "single" if eps < EPS_SAT else "two-phase"
        f.write(f"{eps:>7.3f} {regime:>12} {len(df):>6} "
                f"{fe_peak:>13.2f}% {V_p:>+8.3f} {j_p:>16.2f}\n")
    f.write("\nBloomquist target (CEJ 2026): FE_ADPN = 73–80% at ε_org ≈ 0.15, "
            "j = 70–300 mA/cm²\n")
print(f"wrote {out_txt}")

plt.close("all")
print("\nStage 2 plots complete →", PLOT_DIR)
