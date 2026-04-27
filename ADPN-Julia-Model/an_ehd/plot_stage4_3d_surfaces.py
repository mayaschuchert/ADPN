"""
plot_stage4_3d_surfaces.py — Bloomquist Fig. 5 recreation with v6 model overlay.

Inputs:
  output/stage4/data/stage4_diagnostic.csv   (per-row model FE_ADN, V_cell, etc.)

Outputs:
  output/stage4/plots/stage4_3d_FE_ADN_obs.png       Bloomquist data: We_aq × We_org × j → FE_ADN
  output/stage4/plots/stage4_3d_FE_ADN_model.png     v6 model: We_aq × We_org × j → FE_ADN
  output/stage4/plots/stage4_3d_FE_ADN_resid.png     residual (model − obs)
  output/stage4/plots/stage4_3d_FE_ADN_obs_vs_model.png  side-by-side comparison

Each figure has 3 panels (one per gap). Per v6 §21 panel (o) — the regime-map
visual test. Goal: does the model put the high-FE_ADN region in the same
We_aq / We_org / j corner as the experiment, *independent of magnitude*?

No GPR / interpolation — just colored 3D scatter of the experimental (or model)
points. Bloomquist's Fig. 5 used GPR contour interpolation; that requires
scipy or scikit-learn and is deferred to a v7 plot if reviewers ask for it.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_FN  = os.path.join(HERE, "output", "stage4", "data", "stage4_diagnostic.csv")
EXP_FN   = os.path.join(HERE, "Experimental_data", "bloomquist_data.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_ORDER = (0.25, 0.5, 1.0)


def load_with_weber():
    """Load diagnostic.csv and merge in We_aq / We_org from bloomquist_data.csv."""
    diag = pd.read_csv(DATA_FN)
    exp  = pd.read_csv(EXP_FN)
    keys = ["table", "gap_mm", "Q_total_mL_min", "j_mA_cm2", "phi_AN"]
    we   = exp[keys + ["We_aq", "We_org"]].copy()
    merged = diag.merge(we, on=keys, how="left")
    return merged


def scatter_panel(ax, df, value_col, vmin, vmax, cmap, title):
    """3D scatter at log10(We_aq), log10(We_org), j_mA_cm2 colored by value_col."""
    sub = df[np.isfinite(df[value_col]) & (df["We_aq"] > 0) & (df["We_org"] > 0)]
    sc = ax.scatter(np.log10(sub["We_aq"]),
                    np.log10(sub["We_org"]),
                    sub["j_mA_cm2"],
                    c=sub[value_col].values,
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    s=44, edgecolors="0.3", linewidths=0.4)
    ax.set_xlabel("log₁₀ We_aq")
    ax.set_ylabel("log₁₀ We_org")
    ax.set_zlabel("j [mA/cm²]")
    ax.set_title(title, fontsize=10)
    ax.view_init(elev=22, azim=-58)
    return sc


def make_three_panel_fig(df, value_col, vmin, vmax, cmap, suptitle, fname,
                         cbar_label):
    fig = plt.figure(figsize=(15, 5))
    sc = None
    for k, gap in enumerate(GAP_ORDER):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        gdf = df[df["gap_mm"] == gap]
        sc = scatter_panel(ax, gdf, value_col, vmin, vmax, cmap,
                           f"gap = {gap} mm  (n = {len(gdf)})")
    cbar_ax = fig.add_axes([0.92, 0.18, 0.013, 0.65])
    fig.colorbar(sc, cax=cbar_ax, label=cbar_label)
    fig.suptitle(suptitle, fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.04, right=0.90, wspace=0.10,
                        top=0.92, bottom=0.06)
    path = os.path.join(PLOT_DIR, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


def make_side_by_side(df):
    """Two rows × 3 gaps: row 1 = Bloomquist FE_ADN, row 2 = model FE_ADN."""
    fig = plt.figure(figsize=(15, 9))
    sc = None
    for row, (col, label) in enumerate(
            [("FE_ADN_obs",   "Bloomquist FE_ADN"),
             ("FE_ADN_model", "v6 model FE_ADN")]):
        for k, gap in enumerate(GAP_ORDER):
            idx = row * 3 + k + 1
            ax = fig.add_subplot(2, 3, idx, projection="3d")
            gdf = df[df["gap_mm"] == gap]
            sc = scatter_panel(ax, gdf, col, vmin=0, vmax=80, cmap="turbo",
                               title=f"{label} | gap = {gap} mm")
    cbar_ax = fig.add_axes([0.92, 0.20, 0.013, 0.60])
    fig.colorbar(sc, cax=cbar_ax, label="FE_ADN [%]")
    fig.suptitle("Stage 4 — Bloomquist (top) vs v6 model (bottom): "
                 "FE_ADN over (We_aq, We_org, j)", fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.04, right=0.90, wspace=0.05, hspace=0.18,
                        top=0.94, bottom=0.05)
    path = os.path.join(PLOT_DIR, "stage4_3d_FE_ADN_obs_vs_model.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    df = load_with_weber()
    sub = df[df["converged"] == True].copy()
    print(f"loaded {len(df)} rows ({len(sub)} converged)")

    # Bloomquist FE_ADN (uses all rows — observations are present even when
    # the model failed to converge)
    make_three_panel_fig(df, "FE_ADN_obs", vmin=0, vmax=80, cmap="turbo",
                         suptitle="Stage 4 — Bloomquist FE_ADN (experiment)",
                         fname="stage4_3d_FE_ADN_obs.png",
                         cbar_label="FE_ADN [%]")

    # Model FE_ADN (only converged rows have a meaningful number)
    make_three_panel_fig(sub, "FE_ADN_model", vmin=0, vmax=80, cmap="turbo",
                         suptitle="Stage 4 — v6 model FE_ADN (fitted kinetics)",
                         fname="stage4_3d_FE_ADN_model.png",
                         cbar_label="FE_ADN [%]")

    # Residual (model − obs)
    make_three_panel_fig(sub, "FE_ADN_resid_pp",
                         vmin=-30, vmax=30, cmap="RdBu_r",
                         suptitle="Stage 4 — FE_ADN residual (model − obs)",
                         fname="stage4_3d_FE_ADN_resid.png",
                         cbar_label="ΔFE_ADN [pp]")

    # Side-by-side comparison
    make_side_by_side(sub)


if __name__ == "__main__":
    main()
