"""
plot_fe_vs_j_sweep.py — FE_ADN, FE_PN, FE_TCH, FE_HER vs total current density
for the 6 phi_AN (eps_org) values from stage2m records.
Overlays experimental scatter (Bloomquist data) for FE_ADN, FE_PN, FE_TCH.
FE_HER panel shows model curves only (no experimental data available).

Input:
  output/data/stage2m_records_eo<X>_d50.csv   (6 files)
  ../Experimental_data/bloomquist_data.csv

Output:
  output/stage4_seq/plots/fe_vs_j_sweep.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

HERE     = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "output", "data")
EXP_CSV  = os.path.join(HERE, "..", "Experimental_data", "bloomquist_data.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4_seq", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# Paul Tol muted 6-color sequence (colorblind-safe)
TOL_MUTED = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677"]

EO_VALS = ["0.020", "0.050", "0.090", "0.150", "0.250", "0.300"]

matplotlib.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       10,
    "axes.linewidth":  1.0,
    "axes.labelsize":  10,
    "legend.fontsize": 8.5,
    "figure.dpi":      150,
})


def load_stage2m():
    frames = {}
    for eo in EO_VALS:
        path = os.path.join(DATA_DIR, f"stage2m_records_eo{eo}_d50.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        df["j_mA_cm2"] = df["j_total"] * 0.1   # A m⁻² → mA cm⁻²
        df["phi_AN"]   = float(eo)
        frames[eo] = df
    return frames


def load_exp():
    df = pd.read_csv(EXP_CSV)
    eo_floats = [float(v) for v in EO_VALS]
    mask = df["phi_AN"].apply(lambda x: any(abs(x - e) < 0.005 for e in eo_floats))
    return df[mask].copy()


def style_ax(ax, xmax=None):
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="minor", length=2.5)
    ax.grid(True, linewidth=0.4, alpha=0.45)
    ax.set_xlim(0, xmax)


def main():
    frames = load_stage2m()
    exp    = load_exp()

    if not frames:
        print("No stage2m CSV files found — run stage2m first.")
        return

    # Clip x-axis to experimental current range (default kinetics can produce
    # unrealistic j at high V before fitting; only the experimental range matters)
    xmax = exp["j_mA_cm2"].max() * 1.15 if len(exp) else 350.0

    # (model_col, exp_col or None, ylabel, ylim, show_exp)
    channels = [
        ("FE_ADPN", "FE_ADN_pct", "FE$_{\\mathrm{ADN}}$  [%]", (0, 100), True),
        ("FE_PN",   "FE_PN_pct",  "FE$_{\\mathrm{PN}}$  [%]",  (0,  30), True),
        ("FE_TCH",  "FE_TCH_pct", "FE$_{\\mathrm{TCH}}$  [%]", (0,  30), True),
        ("FE_HER",  None,         "FE$_{\\mathrm{HER}}$  [%]", (0,  30), False),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.0), sharey=False)

    legend_handles = []

    for idx, (model_col, exp_col, ylabel, ylim, show_exp) in enumerate(channels):
        ax = axes[idx]
        for i, eo in enumerate(EO_VALS):
            if eo not in frames:
                continue
            df    = frames[eo]
            col   = TOL_MUTED[i]
            label = f"φ_AN = {float(eo):.3f}"

            if model_col not in df.columns:
                continue

            j_mod  = df["j_mA_cm2"].values
            fe_mod = df[model_col].values
            valid  = np.isfinite(j_mod) & np.isfinite(fe_mod)
            if valid.sum() < 2:
                continue
            order = np.argsort(j_mod[valid])
            ax.plot(j_mod[valid][order], fe_mod[valid][order],
                    color=col, linewidth=1.8, linestyle="solid")

            if show_exp and exp_col is not None:
                sub = exp[np.abs(exp["phi_AN"] - float(eo)) < 0.005]
                if len(sub):
                    ax.scatter(sub["j_mA_cm2"], sub[exp_col],
                               color=col, marker="o", s=35,
                               edgecolors=col, facecolors="none",
                               linewidths=1.3, zorder=5)

            if idx == 0:
                legend_handles.append(
                    mlines.Line2D([], [], color=col, linewidth=1.8, label=label))

        ax.set_xlabel("Total current density  [mA cm⁻²]")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        style_ax(ax, xmax=xmax)

    model_line = mlines.Line2D([], [], color="#555555", linewidth=1.8,
                                linestyle="solid", label="Model (stage2m)")
    exp_dot    = mlines.Line2D([], [], color="#555555", marker="o",
                                markersize=7, linestyle="none",
                                markerfacecolor="none", markeredgewidth=1.3,
                                label="Experiment (Bloomquist)")

    all_handles = legend_handles + [model_line, exp_dot]
    fig.legend(all_handles, [h.get_label() for h in all_handles],
               loc="upper center", ncol=5, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.01), frameon=False)

    fig.suptitle(
        "Faradaic efficiency vs current density — φ_AN sweep (stage2m)",
        y=1.06, fontsize=12,
    )
    fig.tight_layout()

    out = os.path.join(PLOT_DIR, "fe_vs_j_sweep.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
