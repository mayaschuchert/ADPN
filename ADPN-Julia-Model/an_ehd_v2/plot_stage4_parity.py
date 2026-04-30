"""
plot_stage4_parity.py — model vs experiment parity for the v7 Stage 4v3 fit.

Inputs (output/stage4v3/data/):
  stage4a_core_residuals.csv
  stage4b_extended_residuals.csv
  stage4b_holdout_residuals.csv

Outputs (output/stage4v3/plots/):
  stage4v3_parity_FE_ADN.png
  stage4v3_parity_FE_PN.png
  stage4v3_parity_FE_TCH.png
  stage4v3_parity_combined.png   3-panel composite (ADN, PN, TCH)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "output", "stage4v3", "data")
PLOT_DIR = os.path.join(HERE, "output", "stage4v3", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_COLOR  = {0.25: "#d62728", 0.5: "#1f77b4", 1.0: "#2ca02c"}
GAP_MARKER = {0.25: "o",       0.5: "s",       1.0: "^"}

SUBSET_LABEL = {
    "Core":     "Core (training)",
    "Extended": "Extended (high-j)",
    "Holdout":  "Holdout (gap = 0.25 mm)",
}


def load_data():
    frames = []
    for fname, tag in [
        ("stage4a_core_residuals.csv",     "Core"),
        ("stage4b_extended_residuals.csv", "Extended"),
        ("stage4b_holdout_residuals.csv",  "Holdout"),
    ]:
        path = os.path.join(DATA_DIR, fname)
        df = pd.read_csv(path)
        df["subset"] = tag
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def rmse(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return np.sqrt(np.mean(x**2)) if x.size else np.nan


def parity_panel(ax, df, obs_col, model_col, label, vmin, vmax, *,
                 title=None, unit=""):
    sub = df[np.isfinite(df[obs_col]) & np.isfinite(df[model_col])]
    for gap, gap_df in sub.groupby("gap_mm"):
        ax.scatter(gap_df[obs_col], gap_df[model_col],
                   c=GAP_COLOR.get(gap, "k"),
                   marker=GAP_MARKER.get(gap, "o"),
                   s=28, alpha=0.75, edgecolors="none",
                   label=f"gap = {gap} mm")
    ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1, alpha=0.5, label="slope 1")
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect("equal")
    ax.set_xlabel(f"Bloomquist {label}{unit}")
    ax.set_ylabel(f"Model {label}{unit}")
    if title is not None:
        ax.set_title(title)
    lines = []
    for tag in ("Core", "Extended", "Holdout"):
        s = sub[sub["subset"] == tag]
        if len(s) == 0:
            continue
        r = rmse(s[model_col] - s[obs_col])
        lines.append(f"{SUBSET_LABEL[tag]}: RMSE = {r:.2f}{unit}  (n={len(s)})")
    ax.text(0.04, 0.96, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="top", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9, pad=4))
    ax.grid(alpha=0.3)


def main():
    df = load_data()
    print(f"loaded {len(df)} rows")
    print(df["subset"].value_counts())

    channels = [
        ("FE_ADN_obs", "FE_ADN_model", "stage4v3_parity_FE_ADN.png",
         "FE_ADN", 0, 100, " [%]", "(a) FE_ADN"),
        ("FE_PN_obs",  "FE_PN_model",  "stage4v3_parity_FE_PN.png",
         "FE_PN",  0,  50, " [%]", "(b) FE_PN"),
        ("FE_TCH_obs", "FE_TCH_model", "stage4v3_parity_FE_TCH.png",
         "FE_TCH", 0,  50, " [%]", "(c) FE_TCH"),
    ]

    for obs_col, model_col, fname, label, vmin, vmax, unit, panel_title in channels:
        fig, ax = plt.subplots(figsize=(6, 6))
        parity_panel(ax, df, obs_col, model_col, label, vmin, vmax,
                     title=f"{label} — model vs Bloomquist (v7)", unit=unit)
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
        fig.tight_layout()
        path = os.path.join(PLOT_DIR, fname)
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"wrote {path}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    for ax, (obs_col, model_col, _, label, vmin, vmax, unit, panel_title) in \
            zip(axes, channels):
        parity_panel(ax, df, obs_col, model_col, label, vmin, vmax,
                     title=panel_title, unit=unit)
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle("Stage 4v3 — model vs Bloomquist parity (v7 kinetics, TCH added)",
                 y=1.06, fontsize=12)
    fig.tight_layout()
    path = os.path.join(PLOT_DIR, "stage4v3_parity_combined.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
