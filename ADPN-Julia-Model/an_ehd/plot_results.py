"""
plot_results.py — Stage 1 diagnostic plots for the ADPN EHD 1D model.

Generates the Stage 1 subset of the plot list (§19):
  (a) j_r vs V
  (b) FE_r vs V
  (e) c_AN(0) / c_AN,bulk vs j
  (f) φ_l(0) vs j
  + species depth profiles at the most-negative V
  + physicality diagnostics panel
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_DIR  = os.path.join(HERE, "output", "data")
PLOT_DIR  = os.path.join(HERE, "output", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)


# ---------------------------------------------------------------
# 1. Load Stage 1 records
# ---------------------------------------------------------------
rec_files = sorted(glob.glob(os.path.join(DATA_DIR, "stage1_records_eo*.csv")))
if not rec_files:
    raise SystemExit("No Stage 1 records found in " + DATA_DIR)

df = pd.concat([pd.read_csv(f) for f in rec_files], ignore_index=True)
df = df.sort_values("V", ascending=False).reset_index(drop=True)

print(f"Loaded {len(df)} converged points from {len(rec_files)} record files.")
print(f"V range: [{df.V.min():.3f}, {df.V.max():.3f}] V vs SHE")

# j in mA/cm² for plotting (model stores A/m²; 1 A/m² = 0.1 mA/cm²)
df["j_total_mA"] = df["j_total"] * 0.1
df["j1_mA"]      = df["j1"]     * 0.1
df["j2_mA"]      = df["j2"]     * 0.1
df["j3_mA"]      = df["j3"]     * 0.1

# ---------------------------------------------------------------
# Panel A1 (log) + A2 (linear) + B: polarization & FE vs V
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
neg_V = -df.V

ax = axes[0]
ax.semilogy(neg_V, df.j1_mA, "-o", label=r"$j_{\rm ADPN}$", color="C0", ms=3)
ax.semilogy(neg_V, df.j2_mA, "-s", label=r"$j_{\rm PN}$",   color="C1", ms=3)
ax.semilogy(neg_V, df.j3_mA, "-^", label=r"$j_{\rm HER}$",  color="C2", ms=3)
ax.semilogy(neg_V, df.j_total_mA, "k--", label=r"$j_{\rm total}$", lw=1.2)
ax.set_xlabel(r"$-V$ vs SHE  [V]")
ax.set_ylabel(r"$j_r$  [mA cm$^{-2}$]")
ax.set_title("Panel (a1): Polarization — log scale")
ax.legend(loc="lower right", fontsize=9)
ax.grid(True, which="both", ls=":")

ax = axes[1]
ax.plot(neg_V, df.j1_mA, "-o", label=r"$j_{\rm ADPN}$", color="C0", ms=3)
ax.plot(neg_V, df.j2_mA, "-s", label=r"$j_{\rm PN}$",   color="C1", ms=3)
ax.plot(neg_V, df.j3_mA, "-^", label=r"$j_{\rm HER}$",  color="C2", ms=3)
ax.plot(neg_V, df.j_total_mA, "k--", label=r"$j_{\rm total}$", lw=1.2)
ax.set_xlabel(r"$-V$ vs SHE  [V]")
ax.set_ylabel(r"$j_r$  [mA cm$^{-2}$]")
ax.set_title("Panel (a2): Polarization — linear scale")
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, ls=":")

ax = axes[2]
ax.plot(neg_V, df.FE_ADPN, "-o", label=r"FE$_{\rm ADPN}$", color="C0", ms=3)
ax.plot(neg_V, df.FE_PN,   "-s", label=r"FE$_{\rm PN}$",   color="C1", ms=3)
ax.plot(neg_V, df.FE_HER,  "-^", label=r"FE$_{\rm HER}$",  color="C2", ms=3)
ax.set_xlabel(r"$-V$ vs SHE  [V]")
ax.set_ylabel("FE  [%]")
ax.set_title("Panel (b): Faradaic efficiencies")
ax.set_ylim(-2, 102)
ax.legend(loc="center right", fontsize=9)
ax.grid(True, ls=":")

_eps_label = df.eps_org.iloc[0]
_delta_label = df.delta_um.iloc[0]
fig.suptitle(
    rf"Stage 1 — $\varepsilon_{{\rm org}}={_eps_label:.3f}$, "
    rf"$\delta={_delta_label:.0f}\,\mu$m, single-phase reference",
    fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage1_panel_ab_polarization_FE.png"), dpi=160)
print("wrote stage1_panel_ab_polarization_FE.png")

# ---------------------------------------------------------------
# Panel E+F: AN depletion + φ_l(0) vs j_total
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

ax = axes[0]
ax.semilogx(df.j_total_mA.abs() + 1e-8, df.AN_depletion, "-o", color="C3", ms=3)
ax.set_xlabel(r"$|j_{\rm total}|$  [mA cm$^{-2}$]")
ax.set_ylabel(r"$c_{\rm AN}(0) / c_{\rm AN,bulk}$")
ax.set_ylim(-0.05, 1.05)
ax.set_title("Panel (e): AN depletion at electrode")
ax.grid(True, which="both", ls=":")

ax = axes[1]
phi_mV = df.phi_l_surface * 1000.0
ax.semilogx(df.j_total_mA.abs() + 1e-8, phi_mV, "-o", color="C4", ms=3)
ax.set_xlabel(r"$|j_{\rm total}|$  [mA cm$^{-2}$]")
ax.set_ylabel(r"$\phi_\ell(0)$  [mV]")
ax.set_title("Panel (f): electrolyte φ at surface (ohmic penalty)")
ax.grid(True, which="both", ls=":")

fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage1_panel_ef_depletion_phi.png"), dpi=160)
print("wrote stage1_panel_ef_depletion_phi.png")

# ---------------------------------------------------------------
# Species depth profiles at the most-negative V
# ---------------------------------------------------------------
profile_files = sorted(glob.glob(os.path.join(DATA_DIR, "stage1_profile_eo*.csv")))
profiles = []
for pf in profile_files:
    try:
        V_val = float(pf.split("_V")[-1].replace(".csv", ""))
    except ValueError:
        continue
    profiles.append((V_val, pd.read_csv(pf)))
profiles.sort(key=lambda p: -p[0])   # descending V (less negative first)

def plot_one_profile(V_val, dp, fname, title_suffix):
    fig, axes = plt.subplots(2, 3, figsize=(17, 8))
    x_um = dp.x_m * 1e6

    # pH = -log10(c_H / 1000) where c_H is in mol/m³ → mol/L via /1000
    pH = -np.log10(np.maximum(dp.c_H.values / 1000.0, 1e-300))

    ax = axes[0, 0]
    ax.semilogy(x_um, dp.c_H,  label=r"H$^+$")
    ax.semilogy(x_um, dp.c_OH, label=r"OH$^-$")
    ax.set_xlabel(r"$x$ [μm]");  ax.set_ylabel(r"$c$ [mol m$^{-3}$]")
    ax.set_title("H⁺ and OH⁻ profiles")
    ax.legend(); ax.grid(True, which="both", ls=":")

    ax = axes[0, 1]
    ax.plot(x_um, dp.c_H2PO4, label=r"H$_2$PO$_4^-$")
    ax.plot(x_um, dp.c_HPO4,  label=r"HPO$_4^{2-}$")
    ax.plot(x_um, dp.c_PO4,   label=r"PO$_4^{3-}$")
    ax.set_xlabel(r"$x$ [μm]");  ax.set_ylabel(r"$c$ [mol m$^{-3}$]")
    ax.set_title("Phosphate speciation")
    ax.legend(); ax.grid(True, ls=":")

    ax = axes[0, 2]
    ax.plot(x_um, pH, "C6-", lw=1.5)
    ax.axhline(13.0156, color="gray", ls="--", lw=0.8, label="bulk pH = 13.02")
    ax.set_xlabel(r"$x$ [μm]");  ax.set_ylabel("pH")
    ax.set_title(r"pH profile  ($\mathrm{pH} = -\log_{10}(c_{\mathrm{H}^+}\,/\,1000)$)")
    ax.legend(loc="best", fontsize=9); ax.grid(True, ls=":")

    ax = axes[1, 0]
    ax.plot(x_um, dp.c_AN,   label=r"AN",   color="C3")
    ax.plot(x_um, dp.c_ADPN, label=r"ADPN", color="C0")
    ax.plot(x_um, dp.c_PN,   label=r"PN",   color="C1")
    ax.set_xlabel(r"$x$ [μm]");  ax.set_ylabel(r"$c$ [mol m$^{-3}$]")
    ax.set_title("AN / ADPN / PN profiles")
    ax.legend(); ax.grid(True, ls=":")

    ax = axes[1, 1]
    ax.plot(x_um, dp.phi_l * 1000, "k-")
    ax.set_xlabel(r"$x$ [μm]");  ax.set_ylabel(r"$\phi_\ell$ [mV]")
    ax.set_title("Electrolyte potential")
    ax.grid(True, ls=":")

    # Reserve panel (1, 2) — leave blank for now so the layout stays 2×3.
    axes[1, 2].axis("off")

    fig.suptitle(
        f"Stage 1 depth profiles at V = {V_val:+.3f} V  {title_suffix}",
        fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, fname), dpi=160)
    print(f"wrote {fname}")

# Per-V profile plots
for V_val, dp in profiles:
    fname = f"stage1_profile_V{V_val:+.3f}.png".replace("+", "p").replace("-", "m")
    plot_one_profile(V_val, dp, fname, "(single V)")

# Pick a "canonical" profile nearest V=-2.0 V (around the FE peak) for the default plot
if profiles:
    best = min(profiles, key=lambda p: abs(p[0] + 2.0))
    plot_one_profile(best[0], best[1], "stage1_profiles.png",
                     "(nearest to FE peak; in-range)")

# ---------------------------------------------------------------
# Physicality diagnostics
# ---------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
neg_V = -df.V

ax = axes[0, 0]
ax.plot(neg_V, df.pH_bulk, "-o", ms=3)
ax.axhline(13.02, color="gray", ls="--", label="expected bulk pH ≈ 13.02")
ax.set_xlabel(r"$-V$ [V]"); ax.set_ylabel("bulk pH")
ax.set_title("Bulk pH (should stay at 13.02)")
ax.legend(); ax.grid(True, ls=":")

ax = axes[0, 1]
ax.semilogy(neg_V, df.phi_span * 1000, "-o", ms=3)
ax.set_xlabel(r"$-V$ [V]"); ax.set_ylabel(r"$\phi_\ell$ span [mV]")
ax.set_title(r"Electrolyte potential span across DL")
ax.grid(True, which="both", ls=":")

ax = axes[1, 0]
ax.semilogy(neg_V, df.max_eneutr + 1e-20, "-o", ms=3, color="C2")
ax.set_xlabel(r"$-V$ [V]")
ax.set_ylabel(r"max $|\Sigma z_i c_i|$ (9-species)  [mol m$^{-3}$]")
ax.set_title("9-species electroneutrality  (≈ 0 by construction)")
ax.grid(True, which="both", ls=":")

ax = axes[1, 1]
ax.semilogy(neg_V, df.max_Rbuf_surf, "-o", label="surface", ms=3)
ax.semilogy(neg_V, df.max_Rbuf_bulk, "-s", label="bulk",    ms=3)
ax.set_xlabel(r"$-V$ [V]")
ax.set_ylabel(r"max $|R_{\rm buf}|$  [mol m$^{-3}$ s$^{-1}$]")
ax.set_title("Buffer residuals")
ax.legend(); ax.grid(True, which="both", ls=":")

fig.suptitle("Stage 1 physicality diagnostics", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage1_physicality.png"), dpi=160)
print("wrote stage1_physicality.png")

# ---------------------------------------------------------------
# Na⁺ DL accumulation (supplementary)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.semilogx(df.j_total_mA.abs() + 1e-8, df.Na_dev * 100, "-o", ms=3, color="C5")
ax.set_xlabel(r"$|j_{\rm total}|$  [mA cm$^{-2}$]")
ax.set_ylabel(r"max $|c_{\rm Na}(x) - c_{\rm Na,bulk}| / c_{\rm Na,bulk}$  [%]")
ax.set_title(r"Na⁺ reservoir accumulation in the DL (expected; Na⁺ untransported)")
ax.grid(True, which="both", ls=":")
fig.tight_layout()
fig.savefig(os.path.join(PLOT_DIR, "stage1_Na_accumulation.png"), dpi=160)
print("wrote stage1_Na_accumulation.png")

# ---------------------------------------------------------------
# Compact summary table
# ---------------------------------------------------------------
out_txt = os.path.join(PLOT_DIR, "stage1_summary.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    eps_val = df.eps_org.iloc[0]
    delta_val = df.delta_um.iloc[0]
    f.write(f"Stage 1 summary — ε_org = {eps_val:.3f}, δ = {delta_val:.0f} μm, N = 100, stretch = 10\n")
    f.write("=" * 80 + "\n")
    f.write(f"N converged points : {len(df)}\n")
    f.write(f"V range            : [{df.V.min():+.4f}, {df.V.max():+.4f}] V vs SHE\n")
    f.write(f"j_total range      : [{df.j_total.min():.3e}, {df.j_total.max():.3e}] A/m²\n")
    f.write(f"                   = [{df.j_total_mA.min():.3e}, {df.j_total_mA.max():.3e}] mA/cm²\n")
    f.write(f"max FE_ADPN        : {df.FE_ADPN.max():.2f}%  at V = {df.loc[df.FE_ADPN.idxmax(), 'V']:+.4f} V"
            f"   (j = {df.loc[df.FE_ADPN.idxmax(), 'j_total_mA']:.3e} mA/cm²)\n")
    f.write(f"max FE_PN          : {df.FE_PN.max():.2f}%   at V = {df.loc[df.FE_PN.idxmax(), 'V']:+.4f} V\n")
    f.write(f"max FE_HER         : {df.FE_HER.max():.2f}%  at V = {df.loc[df.FE_HER.idxmax(), 'V']:+.4f} V\n")
    f.write(f"pH_bulk range      : [{df.pH_bulk.min():.4f}, {df.pH_bulk.max():.4f}]  (expected 13.0156)\n")
    f.write(f"max |Σ z_i c_i|    : {df.max_eneutr.max():.3e} mol/m³  (9-species electroneutrality)\n")
    f.write(f"max Na DL dev      : {df.Na_dev.max()*100:.3e}%  (physical; Na⁺ not transported)\n")
    f.write(f"min AN_depletion   : {df.AN_depletion.min():.4f}  (c_AN(0)/c_AN,bulk)\n")
    f.write(f"max φ_l span       : {df.phi_span.max()*1000:.3f} mV across DL\n")
    f.write(f"max |R_buf| surf.  : {df.max_Rbuf_surf.max():.3e} mol/m³/s\n")
    f.write(f"max |R_buf| bulk   : {df.max_Rbuf_bulk.max():.3e} mol/m³/s\n")
print(f"wrote {out_txt}")

plt.close("all")
print("All Stage 1 plots written to:", PLOT_DIR)
