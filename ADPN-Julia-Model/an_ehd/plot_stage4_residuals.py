"""
plot_stage4_residuals.py — residual diagnostic plots for the v6 Stage 4 fit.

Inputs:
  output/stage4/data/stage4_diagnostic.csv

Outputs:
  output/stage4/plots/stage4_resid_FE_ADN_vs_j.png       FE_ADN residual vs j, faceted by gap
  output/stage4/plots/stage4_resid_FE_ADN_vs_eps.png     FE_ADN residual vs ε_org, faceted by gap
  output/stage4/plots/stage4_resid_FE_PN_vs_j.png        FE_PN residual vs j, faceted by gap
  output/stage4/plots/stage4_resid_FE_PN_vs_eps.png      FE_PN residual vs ε_org, faceted by gap
  output/stage4/plots/stage4_resid_V_cell_vs_j.png       V_cell residual vs j, faceted by gap

Per v6 §21 panels (k) and (l). The diagnostic value is in the *shape* of the
residual cloud — random scatter ⇒ kinetics fit OK; systematic correlation
⇒ structural model error in that axis.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_FN  = os.path.join(HERE, "output", "stage4", "data", "stage4_diagnostic.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_ORDER  = (0.25, 0.5, 1.0)
SUBSET_COLOR = {
    "Core":          "#1f77b4",
    "Extended-only": "#ff7f0e",
    "Holdout":       "#d62728",
    "Excluded":      "0.6",
}
SUBSET_MARKER = {
    "Core":          "o",
    "Extended-only": "s",
    "Holdout":       "^",
    "Excluded":      "x",
}


def faceted_residual(df, x_col, resid_col, x_label, resid_label, fname, *,
                     ymax=None, title_extra=""):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    sub = df[df["converged"] == True].copy()
    sub = sub[np.isfinite(sub[x_col]) & np.isfinite(sub[resid_col])]

    for ax, gap in zip(axes, GAP_ORDER):
        gdf = sub[sub["gap_mm"] == gap]
        for tag, sdf in gdf.groupby("subset"):
            ax.scatter(sdf[x_col], sdf[resid_col],
                       c=SUBSET_COLOR.get(tag, "k"),
                       marker=SUBSET_MARKER.get(tag, "o"),
                       s=30, alpha=0.8, edgecolors="none",
                       label=tag if ax is axes[0] else None)
        ax.axhline(0, color="k", lw=1, alpha=0.5)
        ax.set_title(f"gap = {gap} mm  (n = {len(gdf)})")
        ax.set_xlabel(x_label)
        if ax is axes[0]:
            ax.set_ylabel(resid_label)
        # RMSE annotation (Core / Extended-only / Holdout)
        annot_lines = []
        for tag in ("Core", "Extended-only", "Holdout"):
            sdf = gdf[gdf["subset"] == tag]
            if len(sdf) == 0:
                continue
            rmse = np.sqrt(np.mean(sdf[resid_col].values ** 2))
            annot_lines.append(f"{tag}: RMSE = {rmse:.2f}  (n={len(sdf)})")
        if annot_lines:
            ax.text(0.04, 0.96, "\n".join(annot_lines),
                    transform=ax.transAxes, ha="left", va="top", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="0.8",
                              alpha=0.9, pad=4))
        if ymax is not None:
            ax.set_ylim(-ymax, ymax)
        ax.grid(alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center",
               ncol=4, fontsize=9, bbox_to_anchor=(0.5, 1.02),
               frameon=False)
    fig.suptitle(f"Stage 4 — {resid_label} vs {x_label}{title_extra}",
                 y=1.06, fontsize=12)
    fig.tight_layout()
    path = os.path.join(PLOT_DIR, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def main():
    df = pd.read_csv(DATA_FN)
    print(f"loaded {len(df)} rows")

    # FE_ADN residuals
    faceted_residual(df, "j_mA_cm2",     "FE_ADN_resid_pp",
                     "j [mA/cm²]",         "FE_ADN residual [pp]",
                     "stage4_resid_FE_ADN_vs_j.png", ymax=50)
    faceted_residual(df, "phi_AN",       "FE_ADN_resid_pp",
                     "ε_org (= φ_AN)",     "FE_ADN residual [pp]",
                     "stage4_resid_FE_ADN_vs_eps.png", ymax=50)

    # FE_PN residuals
    faceted_residual(df, "j_mA_cm2",     "FE_PN_resid_pp",
                     "j [mA/cm²]",         "FE_PN residual [pp]",
                     "stage4_resid_FE_PN_vs_j.png", ymax=20)
    faceted_residual(df, "phi_AN",       "FE_PN_resid_pp",
                     "ε_org (= φ_AN)",     "FE_PN residual [pp]",
                     "stage4_resid_FE_PN_vs_eps.png", ymax=20)

    # V_cell residuals (model − obs, in V)
    df["V_cell_resid_V"] = df["V_cell_pred_V"] - df["V_cell_obs_V"]
    faceted_residual(df, "j_mA_cm2",     "V_cell_resid_V",
                     "j [mA/cm²]",         "V_cell residual [V]",
                     "stage4_resid_V_cell_vs_j.png", ymax=2)


if __name__ == "__main__":
    main()
