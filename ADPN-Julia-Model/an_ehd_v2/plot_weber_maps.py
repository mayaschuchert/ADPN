"""
plot_weber_maps.py — 3D GPR cube maps: We_aq × We_org × j → FE_ADN.
Matches Bloomquist Fig. 5 style: three solid-colored cube faces (GPR surface
on right wall, back wall, and floor) colored by predicted FE_ADN [0–100 %].

Input:  ../Experimental_data/bloomquist_data.csv
Output: output/stage4_seq/plots/weber_maps.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
from sklearn.preprocessing import StandardScaler

HERE     = os.path.dirname(os.path.abspath(__file__))
EXP_CSV  = os.path.join(HERE, "..", "Experimental_data", "bloomquist_data.csv")
PLOT_DIR = os.path.join(HERE, "output", "stage4_seq", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

GAP_ORDER  = [0.25, 0.5, 1.0]
GAP_LABELS = ["0.25 mm", "0.5 mm", "1.0 mm"]

N_GRID = 35          # grid resolution per axis on each face
ELEV, AZIM = 25, 135   # upper-left-behind → left wall (y_max), right wall (x_min), top (z_max)

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size":   14,
    "figure.dpi":  150,
})


# ── GPR fit ──────────────────────────────────────────────────────────────────
def fit_gpr(x, y, z, fe):
    X = np.column_stack([x, y, z])
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    kernel = (
        ConstantKernel(50.0, (1e-1, 1e3))
        * RBF([1.0, 1.0, 1.0], (1e-2, 10.0))
        + WhiteKernel(4.0, (1e-2, 100.0))
    )
    gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5,
                                   normalize_y=True, random_state=0)
    gpr.fit(Xs, fe)
    fe_hat = np.clip(gpr.predict(Xs), 0.0, 100.0)
    ss_res = np.sum((fe - fe_hat) ** 2)
    ss_tot = np.sum((fe - fe.mean()) ** 2)
    r2   = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(ss_res / len(fe)))
    return gpr, sc, rmse, r2


def predict_face(gpr, sc, a_vals, b_vals, fixed_idx, fixed_val):
    """Predict FE_ADN on a 2-D face of the bounding box, clipped to [0, 100]."""
    A, B = np.meshgrid(a_vals, b_vals, indexing="ij")   # (Na, Nb)
    pts = np.zeros((A.size, 3))
    pts[:, fixed_idx] = fixed_val
    other = [i for i in range(3) if i != fixed_idx]
    pts[:, other[0]] = A.ravel()
    pts[:, other[1]] = B.ravel()
    Xs = sc.transform(pts)
    fe = np.clip(gpr.predict(Xs), 0.0, 100.0).reshape(A.shape)
    return A, B, fe


def face_colors(fe_arr, cmap, norm):
    """Per-facet RGBA: average the 4 corner FE values, shape (Na-1, Nb-1, 4)."""
    fe_f = 0.25 * (fe_arr[:-1, :-1] + fe_arr[1:, :-1]
                   + fe_arr[:-1, 1:] + fe_arr[1:, 1:])
    return cmap(norm(fe_f))


def decade_ticks(vals):
    lo, hi = int(np.floor(vals.min())), int(np.ceil(vals.max()))
    return list(range(lo, hi + 1))


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    df = pd.read_csv(EXP_CSV)
    df = df[df["phi_AN"] > 0.01].copy()
    for col in ("We_aq", "We_org", "j_mA_cm2", "FE_ADN_pct"):
        df = df[np.isfinite(df[col]) & (df[col] > 0)].copy()
    df["lwaq"]  = np.log10(df["We_aq"])
    df["lworg"] = np.log10(df["We_org"])

    cmap = plt.cm.turbo
    norm = Normalize(vmin=0, vmax=80)

    fig = plt.figure(figsize=(21, 8))

    for k, (gap, glabel) in enumerate(zip(GAP_ORDER, GAP_LABELS)):
        gdf = df[np.isclose(df["gap_mm"], gap, atol=0.01)].copy().reset_index(drop=True)
        ax  = fig.add_subplot(1, 3, k + 1, projection="3d")

        x  = gdf["lwaq"].values       # log10(We_aq)
        y  = gdf["lworg"].values      # log10(We_org)
        z  = gdf["j_mA_cm2"].values
        fe = gdf["FE_ADN_pct"].values

        gpr, sc, rmse, r2 = fit_gpr(x, y, z, fe)
        print(f"  gap={gap}: n={len(gdf)}  RMSE={rmse:.2f}  R²={r2:.3f}")

        waq_g  = np.linspace(x.min(), x.max(), N_GRID)
        worg_g = np.linspace(y.min(), y.max(), N_GRID)
        j_g    = np.linspace(z.min(), z.max(), N_GRID)

        # Features: [log_waq=0, log_worg=1, j=2]
        # Camera at azim=135, elev=25 → negative x, positive y, positive z.
        # Three outer cube faces visible (matching Bloomquist Fig. 5):
        #   • Left wall  (y = y_max, We_org at max) — image-left,  grid over (waq × j)
        #   • Right wall (x = x_min, We_aq at min)  — image-right, grid over (worg × j)
        #   • Top        (z = z_max, j at max)       — grid over (waq × worg)

        # ── Left wall (image-left): y = y_max ─────────────────────────────
        WQ, JJ, FE1 = predict_face(gpr, sc, waq_g, j_g, 1, y.max())
        YY1 = np.full_like(WQ, y.max())
        ax.plot_surface(WQ, YY1, JJ,
                        facecolors=face_colors(FE1, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)

        # ── Right wall (image-right): x = x_min ───────────────────────────
        WO, JJ, FE2 = predict_face(gpr, sc, worg_g, j_g, 0, x.min())
        XX2 = np.full_like(WO, x.min())
        ax.plot_surface(XX2, WO, JJ,
                        facecolors=face_colors(FE2, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)

        # ── Top: z = z_max ────────────────────────────────────────────────
        WQ, WO, FE3 = predict_face(gpr, sc, waq_g, worg_g, 2, z.max())
        ZZ3 = np.full_like(WQ, z.max())
        ax.plot_surface(WQ, WO, ZZ3,
                        facecolors=face_colors(FE3, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)

        # ── Raw data points overlaid ───────────────────────────────────────
        ax.scatter(x, y, z, c=fe, cmap=cmap, norm=norm,
                   s=28, edgecolors="k", linewidths=0.4, zorder=5)

        # ── Axes ───────────────────────────────────────────────────────────
        title = glabel
        if np.isfinite(rmse):
            title += f"\nRMSE = {rmse:.2f}   R² = {r2:.3f}"
        ax.set_title(title, fontsize=15, pad=8)
        ax.set_xlabel("We$_\\mathrm{aqueous}$",  fontsize=14, labelpad=14)
        ax.set_ylabel("We$_\\mathrm{organic}$",  fontsize=14, labelpad=14)
        ax.set_zlabel("Current Density (mA cm$^{-2}$)", fontsize=14, labelpad=14)

        xticks = decade_ticks(x)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"$10^{{{t}}}$" for t in xticks], fontsize=12)
        yticks = decade_ticks(y)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"$10^{{{t}}}$" for t in yticks], fontsize=12)
        ax.tick_params(axis="z", labelsize=12)

        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_zlim(z.min(), z.max())
        ax.set_box_aspect([1, 1, 1])   # force cubic panel
        ax.view_init(elev=ELEV, azim=AZIM)

    # ── Colorbar ──────────────────────────────────────────────────────────────
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.92, 0.18, 0.012, 0.64])
    cb = fig.colorbar(sm, cax=cbar_ax, label="FE$_{\\mathrm{ADN}}$  [%]")
    cb.ax.tick_params(labelsize=12)
    cb.set_label("FE$_{\\mathrm{ADN}}$  [%]", fontsize=14)

    fig.suptitle(
        "Weber-number regime map  (Bloomquist data)\n"
        "GPR surface: We$_{\\mathrm{aq}}$ × We$_{\\mathrm{org}}$ × j  →  "
        "FE$_{\\mathrm{ADN}}$ [0–100 %]",
        fontsize=15, y=1.01,
    )
    fig.subplots_adjust(left=0.04, right=0.90, wspace=0.10, top=0.88, bottom=0.06)

    out = os.path.join(PLOT_DIR, "weber_maps.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
