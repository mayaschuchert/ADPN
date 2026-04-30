"""
plot_stage4seq_residuals.py — residual diagnostic plots for the v7 stage4_seq fit.

Inputs:
  output/stage4_seq/data/stage4seq_parity.csv

Outputs (6 plots in output/stage4_seq/plots/):
  resid_FE_ADN_vs_j.png        FE_ADN residual vs j,      faceted by gap
  resid_FE_ADN_vs_eps.png      FE_ADN residual vs ε_org,  faceted by gap
  resid_FE_PN_vs_j.png         FE_PN  residual vs j,      faceted by gap
  resid_FE_PN_vs_eps.png       FE_PN  residual vs ε_org,  faceted by gap
  resid_FE_TCH_vs_j.png        FE_TCH residual vs j,      faceted by gap
  resid_FE_TCH_vs_eps.png      FE_TCH residual vs ε_org,  faceted by gap

Random scatter ⇒ fit OK; systematic trend ⇒ structural model error.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_FN  = os.path.join(HERE, "output", "stage4_seq", "data", "stage4seq_parity.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4_seq", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_ORDER = [0.25, 0.5, 1.0]

SET_COLOR = {
    "core":     "#004488",
    "extended": "#1166BB",
    "holdout":  "#008877",
    "other":    "0.6",
}
SET_MARKER = {
    "core":     "o",
    "extended": "s",
    "holdout":  "^",
    "other":    "x",
}
SET_LABEL = {
    "core":     "Core",
    "extended": "Extended",
    "holdout":  "Holdout",
    "other":    "Other",
}


def faceted_residual(df, x_col, resid_col, x_label, resid_label, fname, ymax=None):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)

    sub = df[df["converged"] == 1].copy()
    sub = sub[np.isfinite(sub[x_col]) & np.isfinite(sub[resid_col])]

    for ax, gap in zip(axes, GAP_ORDER):
        gdf = sub[np.isclose(sub["gap_mm"], gap, atol=0.01)]

        for tag, sdf in gdf.groupby("set"):
            ax.scatter(
                sdf[x_col], sdf[resid_col],
                c=SET_COLOR.get(tag, "k"),
                marker=SET_MARKER.get(tag, "o"),
                s=30, alpha=0.85, edgecolors="none",
                label=SET_LABEL.get(tag, tag) if ax is axes[0] else None,
            )

        ax.axhline(0, color="k", lw=1.0, alpha=0.5)
        ax.set_title(f"gap = {gap} mm  (n = {len(gdf)})")
        ax.set_xlabel(x_label)
        if ax is axes[0]:
            ax.set_ylabel(resid_label)
        ax.grid(alpha=0.30)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

        # RMSE annotation per subset
        annot = []
        for tag in ("core", "extended", "holdout"):
            sdf = gdf[gdf["set"] == tag]
            if len(sdf) == 0:
                continue
            rmse = np.sqrt(np.mean(sdf[resid_col].values ** 2))
            annot.append(f"{SET_LABEL[tag]}: RMSE = {rmse:.2f}  (n={len(sdf)})")
        if annot:
            ax.text(
                0.04, 0.96, "\n".join(annot),
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.90, pad=4),
            )

        if ymax is not None:
            ax.set_ylim(-ymax, ymax)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle(
        f"Stage 4 seq — {resid_label} vs {x_label}",
        y=1.06, fontsize=12,
    )
    fig.tight_layout()

    path = os.path.join(PLOT_DIR, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def main():
    if not os.path.exists(DATA_FN):
        print(f"Missing: {DATA_FN}")
        print("Run gen_stage4seq_parity.jl first.")
        return

    df = pd.read_csv(DATA_FN)
    print(f"Loaded {len(df)} rows  ({(df['converged']==1).sum()} converged)")

    df["FE_ADN_resid"] = df["FE_ADN_pred"] - df["FE_ADN_obs"]
    df["FE_PN_resid"]  = df["FE_PN_pred"]  - df["FE_PN_obs"]
    df["FE_TCH_resid"] = df["FE_TCH_pred"] - df["FE_TCH_obs"]

    faceted_residual(df, "j_mA_cm2", "FE_ADN_resid",
                     "j  [mA cm⁻²]",   "FE_ADN residual  [pp]",
                     "resid_FE_ADN_vs_j.png",   ymax=50)
    faceted_residual(df, "phi_AN",   "FE_ADN_resid",
                     "ε_org  (φ_AN)",   "FE_ADN residual  [pp]",
                     "resid_FE_ADN_vs_eps.png",  ymax=50)

    faceted_residual(df, "j_mA_cm2", "FE_PN_resid",
                     "j  [mA cm⁻²]",   "FE_PN residual  [pp]",
                     "resid_FE_PN_vs_j.png",    ymax=20)
    faceted_residual(df, "phi_AN",   "FE_PN_resid",
                     "ε_org  (φ_AN)",   "FE_PN residual  [pp]",
                     "resid_FE_PN_vs_eps.png",   ymax=20)

    faceted_residual(df, "j_mA_cm2", "FE_TCH_resid",
                     "j  [mA cm⁻²]",   "FE_TCH residual  [pp]",
                     "resid_FE_TCH_vs_j.png",   ymax=20)
    faceted_residual(df, "phi_AN",   "FE_TCH_resid",
                     "ε_org  (φ_AN)",   "FE_TCH residual  [pp]",
                     "resid_FE_TCH_vs_eps.png",  ymax=20)


if __name__ == "__main__":
    main()
