"""
plot_species_profiles.py — boundary-layer species profiles, blue colorblind-friendly.
2×3 subplot layout:  AN | Products | phi_l
                     Phosphate buffer | H⁺/OH⁻ | (empty → overall legend)
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HERE    = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "output", "data")
FIG_DIR = os.path.join(HERE, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

OP_FILE   = os.path.join(OUT_DIR, "profile_operating.csv")
FB_FILES  = sorted(f for f in
                   [os.path.join(OUT_DIR, f"stage2m_profile_eo0.{e}_d50_V-1.520.csv")
                    for e in ("020","050","090","150","250","300")]
                   if os.path.exists(f))

if os.path.exists(OP_FILE):
    profile_file = OP_FILE
    subtitle = "Operating point  (gap=0.5 mm, Q=2 mL/min, φ_AN=0.29, j=800 A m⁻²)"
else:
    profile_file = FB_FILES[-1] if FB_FILES else None
    subtitle = "Stage 2m equilibrium profile"

if profile_file is None:
    raise FileNotFoundError("No profile CSV found.")

print(f"Plotting: {os.path.basename(profile_file)}")
data = np.genfromtxt(profile_file, delimiter=",", names=True)

x_m  = data["x_m"]
x_um = x_m * 1e6

def safe_col(name):
    """Return column; replace inf/nan with NaN."""
    v = data[name].copy().astype(float)
    v[~np.isfinite(v)] = np.nan
    return v

c_H     = safe_col("c_H")
c_OH    = safe_col("c_OH")
c_H2PO4 = safe_col("c_H2PO4")
c_HPO4  = safe_col("c_HPO4")
c_PO4   = safe_col("c_PO4")
c_AN    = safe_col("c_AN")
c_ADPN  = safe_col("c_ADPN")
c_PN    = safe_col("c_PN")
c_TCH   = safe_col("c_TCH")
phi_l   = safe_col("phi_l")

def norm(v):
    bv = v[np.isfinite(v)][-1]
    return v / bv if bv != 0 else v

# ── colorblind-safe blues (Paul Tol muted/bright) ──────────────────────────
B1 = "#004488"   # deep blue
B2 = "#1166BB"   # cobalt
B3 = "#3388DD"   # royal blue
B4 = "#55AAEE"   # sky blue
B5 = "#77CCFF"   # powder blue
T1 = "#008877"   # teal (blue-green, still accessible)
T2 = "#00BBAA"   # cyan-teal

matplotlib.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      10,
    "axes.linewidth": 1.0,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "figure.dpi":     150,
})

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle(f"Boundary-layer species profiles\n{subtitle}", fontsize=12)

def style_ax(ax, xlabel=True):
    ax.set_xlim(left=0, right=x_um.max())
    if xlabel:
        ax.set_xlabel("Distance from electrode surface  (µm)")
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="minor", length=2.5)
    ax.grid(True, linewidth=0.4, alpha=0.5)

# ── (0,0) AN — main reactant ───────────────────────────────────────────────
ax = axes[0, 0]
ax.plot(x_um, norm(c_AN), color=B1, linewidth=2.2, label="AN")
ax.axhline(1.0, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)
ax.set_ylabel("c / c_bulk")
ax.set_title("Adenosine (AN) — main reactant")
ax.legend(loc="lower right")
style_ax(ax, xlabel=False)

# ── (0,1) Products — ADPN, PN, TCH ────────────────────────────────────────
ax = axes[0, 1]
for (vec, label, color, ls) in [
        (c_ADPN, "ADPN", B2, "solid"),
        (c_PN,   "PN",   B4, "dashed"),
        (c_TCH,  "TCH",  T1, "dotted"),
]:
    ax.plot(x_um, vec, color=color, linewidth=1.8, linestyle=ls, label=label)
ax.set_ylabel("Concentration  (mol m⁻³)")
ax.set_title("Electrochemical products")
ax.legend(loc="upper right")
style_ax(ax, xlabel=False)

# ── (0,2) φ_l potential ───────────────────────────────────────────────────
ax = axes[0, 2]
ax.plot(x_um, phi_l * 1e3, color=B1, linewidth=2.0)
ax.axhline(0, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)
ax.set_ylabel("Liquid potential  φ_l  (mV)")
ax.set_title("Potential profile")
style_ax(ax, xlabel=False)

# ── (1,0) Phosphate buffer ─────────────────────────────────────────────────
ax = axes[1, 0]
for (vec, label, color, ls) in [
        (c_H2PO4, "H₂PO₄⁻", B2, "solid"),
        (c_HPO4,  "HPO₄²⁻",  B3, "dashed"),
        (c_PO4,   "PO₄³⁻",   B5, "dotted"),
]:
    ax.plot(x_um, norm(vec), color=color, linewidth=1.8, linestyle=ls, label=label)
ax.axhline(1.0, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)
ax.set_ylabel("c / c_bulk")
ax.set_title("Phosphate buffer species")
ax.legend(loc="lower right")
style_ax(ax)

# ── (1,1) H⁺ / OH⁻ ────────────────────────────────────────────────────────
ax = axes[1, 1]
for (vec, label, color, ls) in [
        (c_H,  "H⁺",  B2, "solid"),
        (c_OH, "OH⁻", T1, "dashed"),
]:
    ax.plot(x_um, norm(vec), color=color, linewidth=1.8, linestyle=ls, label=label)
ax.axhline(1.0, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)
ax.set_ylabel("c / c_bulk")
ax.set_title("Acid-base species")
ax.legend(loc="best")
style_ax(ax)

# ── (1,2) pH profile ───────────────────────────────────────────────────────
ax = axes[1, 2]
pH = 3.0 - np.log10(np.where(c_H > 0, c_H, np.nan))   # c_H in mol/m³, bulk pH≈3
ax.plot(x_um, pH, color=B1, linewidth=2.0)
ax.set_ylabel("pH")
ax.set_title("Local pH profile")
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
style_ax(ax)

fig.tight_layout(rect=[0, 0, 1, 0.94])

out_path = os.path.join(FIG_DIR, "species_profiles.png")
fig.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"Saved: {out_path}")
plt.show()
