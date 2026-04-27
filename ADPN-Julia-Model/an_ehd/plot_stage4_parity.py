"""
plot_stage4_parity.py — model vs experiment parity for the v6 Stage 4 fit.

Inputs:
  output/stage4/data/stage4_diagnostic.csv   (162 rows × 21 cols)

Outputs:
  output/stage4/plots/stage4_parity_FE_ADN.png   FE_ADN model vs obs, colored by gap
  output/stage4/plots/stage4_parity_FE_PN.png    FE_PN  model vs obs, colored by gap
  output/stage4/plots/stage4_parity_V_cell.png   V_cell model vs back-derived from EP/PR/j
  output/stage4/plots/stage4_parity_combined.png 3-panel composite

Reads only diagnostic.csv. No model re-runs. Per v6 §21 i, j and the V_cell
panel added by the §17 cell-voltage decomposition.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_FN  = os.path.join(HERE, "output", "stage4", "data", "stage4_diagnostic.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_COLOR  = {0.25: "#d62728", 0.5: "#1f77b4", 1.0: "#2ca02c"}
GAP_MARKER = {0.25: "o",       0.5: "s",       1.0: "^"}

SUBSET_LABEL = {
    "Core":          "Core (training)",
    "Extended-only": "Extended-only (high-j)",
    "Holdout":       "Holdout (gap = 0.25 mm)",
    "Excluded":      "Excluded (ε_org < 0.04)",
}


def rmse(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return np.sqrt(np.mean(x**2)) if x.size else np.nan


def parity_panel(ax, df, x_col, y_col, label, vmin, vmax, *, only_converged=True,
                 title=None, unit=""):
    sub = df[df["converged"] == True] if only_converged else df
    sub = sub[np.isfinite(sub[x_col]) & np.isfinite(sub[y_col])]
    for gap, gap_df in sub.groupby("gap_mm"):
        ax.scatter(gap_df[x_col], gap_df[y_col],
                   c=GAP_COLOR.get(gap, "k"),
                   marker=GAP_MARKER.get(gap, "o"),
                   s=28, alpha=0.75,
                   edgecolors="none",
                   label=f"gap = {gap} mm")
    ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1, alpha=0.5,
            label="slope 1")
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect("equal")
    ax.set_xlabel(f"Bloomquist {label}{unit}")
    ax.set_ylabel(f"Model {label}{unit}")
    if title is not None:
        ax.set_title(title)
    # RMSE summary in upper-left
    lines = []
    for tag in ("Core", "Extended-only", "Holdout"):
        s = sub[sub["subset"] == tag]
        if len(s) == 0:
            continue
        r = rmse(s[y_col] - s[x_col])
        lines.append(f"{SUBSET_LABEL[tag]}: RMSE = {r:.2f}{unit}  (n={len(s)})")
    ax.text(0.04, 0.96, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="top", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9, pad=4))
    ax.grid(alpha=0.3)


def main():
    df = pd.read_csv(DATA_FN)
    print(f"loaded {len(df)} rows from {DATA_FN}")
    print("subset counts:")
    print(df["subset"].value_counts())

    # --- Single-panel figs ---
    for x_col, y_col, fname, label, vmin, vmax, unit in [
        ("FE_ADN_obs", "FE_ADN_model",
         "stage4_parity_FE_ADN.png", "FE_ADN", 0,  100, " [%]"),
        ("FE_PN_obs",  "FE_PN_model",
         "stage4_parity_FE_PN.png",  "FE_PN",  0,  50,  " [%]"),
        ("V_cell_obs_V", "V_cell_pred_V",
         "stage4_parity_V_cell.png", "V_cell", 1.5, 6.0, " [V]"),
    ]:
        fig, ax = plt.subplots(figsize=(6, 6))
        parity_panel(ax, df, x_col, y_col, label, vmin, vmax,
                     title=f"{label} — model vs Bloomquist", unit=unit)
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
        fig.tight_layout()
        path = os.path.join(PLOT_DIR, fname)
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"wrote {path}")

    # --- 3-panel composite ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    parity_panel(axes[0], df, "FE_ADN_obs",   "FE_ADN_model",
                 "FE_ADN",  0,   100, title="(a) FE_ADN", unit=" [%]")
    parity_panel(axes[1], df, "FE_PN_obs",    "FE_PN_model",
                 "FE_PN",   0,   50,  title="(b) FE_PN",  unit=" [%]")
    parity_panel(axes[2], df, "V_cell_obs_V", "V_cell_pred_V",
                 "V_cell",  1.5, 6.0, title="(c) V_cell", unit=" [V]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle("Stage 4 fit — model vs Bloomquist parity (v6.0 kinetics-only)",
                 y=1.06, fontsize=12)
    fig.tight_layout()
    path = os.path.join(PLOT_DIR, "stage4_parity_combined.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
