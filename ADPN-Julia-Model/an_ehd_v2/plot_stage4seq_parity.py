"""
plot_stage4seq_parity.py — parity plot for Stage 4 sequential fit (v7).

Input:  output/stage4_seq/data/stage4seq_parity.csv
Output: output/stage4_seq/plots/stage4seq_parity_combined.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "output", "stage4_seq", "data")
PLOT_DIR = os.path.join(HERE, "output", "stage4_seq", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

CSV = os.path.join(DATA_DIR, "stage4seq_parity.csv")

# Paul Tol colorblind-safe palette (blue family + teal/orange accent)
SET_STYLE = {
    "core":     dict(color="#004488", marker="o", label="Core (training)"),
    "extended": dict(color="#1166BB", marker="s", label="Extended"),
    "holdout":  dict(color="#008877", marker="^", label="Holdout"),
    "other":    dict(color="#aaaaaa", marker="x", label="Other"),
}

matplotlib.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      10,
    "axes.linewidth": 1.0,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "figure.dpi":     150,
})


def rmse(resid):
    r = np.asarray(resid, dtype=float)
    r = r[np.isfinite(r)]
    return np.sqrt(np.mean(r**2)) if r.size else np.nan


def parity_panel(ax, df, obs_col, pred_col, vmin, vmax, title):
    ok = df[np.isfinite(df[obs_col]) & np.isfinite(df[pred_col])]
    for sname, style in SET_STYLE.items():
        sub = ok[ok["set"] == sname]
        if len(sub) == 0:
            continue
        ax.scatter(sub[obs_col], sub[pred_col],
                   color=style["color"], marker=style["marker"],
                   s=30, alpha=0.80, edgecolors="none",
                   label=style["label"])

    # parity line
    ax.plot([vmin, vmax], [vmin, vmax], color="#555555", lw=1.0,
            linestyle="--", zorder=0)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect("equal")
    ax.set_xlabel(f"Observed {title}  [%]")
    ax.set_ylabel(f"Predicted {title}  [%]")
    ax.set_title(title)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="minor", length=2.5)
    ax.grid(True, linewidth=0.4, alpha=0.45)

    # RMSE annotations
    lines = []
    for sname, style in SET_STYLE.items():
        sub = ok[ok["set"] == sname]
        if len(sub) == 0:
            continue
        r = rmse(sub[pred_col] - sub[obs_col])
        lines.append(f"{style['label']}: RMSE={r:.1f}%  (n={len(sub)})")
    ax.text(0.04, 0.97, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="top", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.90, pad=3))


def main():
    df = pd.read_csv(CSV)
    print(f"Loaded {len(df)} rows from {CSV}")
    print(df["set"].value_counts())

    channels = [
        ("FE_ADN_obs", "FE_ADN_pred", 0, 100, "FE$_{\\mathrm{ADN}}$"),
        ("FE_PN_obs",  "FE_PN_pred",  0,  50, "FE$_{\\mathrm{PN}}$"),
        ("FE_TCH_obs", "FE_TCH_pred", 0,  50, "FE$_{\\mathrm{TCH}}$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, (obs, pred, vmin, vmax, title) in zip(axes, channels):
        parity_panel(ax, df, obs, pred, vmin, vmax, title)

    # shared legend above panels
    handles, labels_ = [], []
    for sname, style in SET_STYLE.items():
        if (df["set"] == sname).any():
            handles.append(
                matplotlib.lines.Line2D([], [], color=style["color"],
                                        marker=style["marker"], linestyle="none",
                                        markersize=7, label=style["label"]))
            labels_.append(style["label"])
    handles.append(matplotlib.lines.Line2D([], [], color="#555555",
                                            linestyle="--", lw=1.0,
                                            label="Parity (slope 1)"))
    labels_.append("Parity (slope 1)")

    fig.legend(handles, labels_, loc="upper center", ncol=len(handles),
               fontsize=9, bbox_to_anchor=(0.5, 1.01), frameon=False)
    fig.suptitle("Stage 4 sequential fit — model vs experiment (v7 kinetics)",
                 y=1.06, fontsize=12)
    fig.tight_layout()

    out = os.path.join(PLOT_DIR, "stage4seq_parity_combined.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    import matplotlib.lines
    main()
