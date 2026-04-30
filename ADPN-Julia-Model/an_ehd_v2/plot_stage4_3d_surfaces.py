"""
plot_stage4_3d_surfaces.py — 3-D scatter regime maps for the v7 Stage 4v3 fit.

Inputs:
  output/stage4v3/data/stage4a_core_residuals.csv
  output/stage4v3/data/stage4b_extended_residuals.csv
  output/stage4v3/data/stage4b_holdout_residuals.csv
  ../Experimental_data/bloomquist_data.csv   (for We_aq, We_org)

Outputs (output/stage4v3/plots/):
  stage4v3_3d_FE_ADN_obs.png
  stage4v3_3d_FE_ADN_model.png
  stage4v3_3d_FE_ADN_resid.png
  stage4v3_3d_FE_ADN_obs_vs_model.png
  stage4v3_3d_FE_TCH_obs.png
  stage4v3_3d_FE_TCH_model.png
  stage4v3_3d_FE_TCH_resid.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "output", "stage4v3", "data")
EXP_FN   = os.path.join(HERE, "..", "Experimental_data", "bloomquist_data.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4v3", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_ORDER = (0.25, 0.5, 1.0)


def load_data():
    frames = []
    for fname, tag in [
        ("stage4a_core_residuals.csv",     "Core"),
        ("stage4b_extended_residuals.csv", "Extended"),
        ("stage4b_holdout_residuals.csv",  "Holdout"),
    ]:
        df = pd.read_csv(os.path.join(DATA_DIR, fname))
        df["subset"] = tag
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)

    exp  = pd.read_csv(EXP_FN)
    keys = ["table", "gap_mm", "Q_total_mL_min", "j_mA_cm2", "phi_AN"]
    we   = exp[keys + ["We_aq", "We_org"]].copy()
    return merged.merge(we, on=keys, how="left")


def scatter_panel(ax, df, value_col, vmin, vmax, cmap, title):
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


def make_side_by_side_adn(df):
    fig = plt.figure(figsize=(15, 9))
    sc = None
    for row, (col, label) in enumerate(
            [("FE_ADN_obs",   "Bloomquist FE_ADN"),
             ("FE_ADN_model", "v7 model FE_ADN")]):
        for k, gap in enumerate(GAP_ORDER):
            idx = row * 3 + k + 1
            ax = fig.add_subplot(2, 3, idx, projection="3d")
            gdf = df[df["gap_mm"] == gap]
            sc = scatter_panel(ax, gdf, col, vmin=0, vmax=80, cmap="turbo",
                               title=f"{label} | gap = {gap} mm")
    cbar_ax = fig.add_axes([0.92, 0.20, 0.013, 0.60])
    fig.colorbar(sc, cax=cbar_ax, label="FE_ADN [%]")
    fig.suptitle("Stage 4v3 — Bloomquist (top) vs v7 model (bottom): "
                 "FE_ADN over (We_aq, We_org, j)", fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.04, right=0.90, wspace=0.05, hspace=0.18,
                        top=0.94, bottom=0.05)
    path = os.path.join(PLOT_DIR, "stage4v3_3d_FE_ADN_obs_vs_model.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    df = load_data()
    print(f"loaded {len(df)} rows")

    make_three_panel_fig(df, "FE_ADN_obs", vmin=0, vmax=80, cmap="turbo",
                         suptitle="Stage 4v3 — Bloomquist FE_ADN (experiment)",
                         fname="stage4v3_3d_FE_ADN_obs.png",
                         cbar_label="FE_ADN [%]")
    make_three_panel_fig(df, "FE_ADN_model", vmin=0, vmax=80, cmap="turbo",
                         suptitle="Stage 4v3 — v7 model FE_ADN (fitted kinetics)",
                         fname="stage4v3_3d_FE_ADN_model.png",
                         cbar_label="FE_ADN [%]")
    make_three_panel_fig(df, "FE_ADN_resid_pp",
                         vmin=-30, vmax=30, cmap="RdBu_r",
                         suptitle="Stage 4v3 — FE_ADN residual (model − obs)",
                         fname="stage4v3_3d_FE_ADN_resid.png",
                         cbar_label="ΔFE_ADN [pp]")
    make_side_by_side_adn(df)

    make_three_panel_fig(df, "FE_TCH_obs", vmin=0, vmax=30, cmap="plasma",
                         suptitle="Stage 4v3 — Bloomquist FE_TCH (experiment)",
                         fname="stage4v3_3d_FE_TCH_obs.png",
                         cbar_label="FE_TCH [%]")
    make_three_panel_fig(df, "FE_TCH_model", vmin=0, vmax=30, cmap="plasma",
                         suptitle="Stage 4v3 — v7 model FE_TCH (fitted kinetics)",
                         fname="stage4v3_3d_FE_TCH_model.png",
                         cbar_label="FE_TCH [%]")
    make_three_panel_fig(df, "FE_TCH_resid_pp",
                         vmin=-15, vmax=15, cmap="RdBu_r",
                         suptitle="Stage 4v3 — FE_TCH residual (model − obs)",
                         fname="stage4v3_3d_FE_TCH_resid.png",
                         cbar_label="ΔFE_TCH [pp]")


if __name__ == "__main__":
    main()
