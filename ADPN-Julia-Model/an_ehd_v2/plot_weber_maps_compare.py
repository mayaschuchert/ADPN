"""
plot_weber_maps_compare.py — generate one Weber-map PNG per (θ candidate)
so they can be visually compared side-by-side.

Reads:
  output/forward_v3/data/parity.csv             → output/forward_v3/data/weber_maps_v3.png
  output/forward_stage4seq_LB/data/parity.csv   → output/forward_stage4seq_LB/data/weber_maps_stage4seq_LB.png

The classproject Run-2 map at:
  ../an_ehd_classproject/output/plots/weber_maps.png
already exists; this script doesn't touch it.
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

HERE = os.path.dirname(os.path.abspath(__file__))

CASES = [
    ("forward_v3",           "v3 theta  (joint LM, post-fix, n_PN=1.64, n_TCH=1.55)",
     "weber_maps_v3.png"),
    ("forward_stage4seq_LB", "stage4_seq LB-pinned theta  (n_PN=0.50@LB, n_TCH=1.00@LB)",
     "weber_maps_stage4seq_LB.png"),
]

GAP_ORDER  = [0.25, 0.5, 1.0]
GAP_LABELS = ["0.25 mm", "0.5 mm", "1.0 mm"]
N_GRID = 35
ELEV, AZIM = 25, 135

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size":   14,
    "figure.dpi":  150,
})


def fit_gpr(x, y, z, fe):
    X = np.column_stack([x, y, z])
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    kernel = (ConstantKernel(50.0, (1e-1, 1e3))
              * RBF([1.0, 1.0, 1.0], (1e-2, 10.0))
              + WhiteKernel(4.0, (1e-2, 100.0)))
    gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5,
                                   normalize_y=True, random_state=0)
    gpr.fit(Xs, fe)
    return gpr, sc


def predict_face(gpr, sc, a_vals, b_vals, fixed_idx, fixed_val):
    A, B = np.meshgrid(a_vals, b_vals, indexing="ij")
    pts = np.zeros((A.size, 3))
    pts[:, fixed_idx] = fixed_val
    other = [i for i in range(3) if i != fixed_idx]
    pts[:, other[0]] = A.ravel()
    pts[:, other[1]] = B.ravel()
    Xs = sc.transform(pts)
    fe = np.clip(gpr.predict(Xs), 0.0, 100.0).reshape(A.shape)
    return A, B, fe


def face_colors(fe_arr, cmap, norm):
    fe_f = 0.25 * (fe_arr[:-1, :-1] + fe_arr[1:, :-1]
                   + fe_arr[:-1, 1:] + fe_arr[1:, 1:])
    return cmap(norm(fe_f))


def decade_ticks(vals):
    lo, hi = int(np.floor(vals.min())), int(np.ceil(vals.max()))
    return list(range(lo, hi + 1))


def make_one_map(parity_csv: str, label: str, out_png: str):
    df = pd.read_csv(parity_csv)
    df = df[df["converged"] == 1].copy()
    df = df[df["phi_AN"] > 0.01].copy()
    for col in ("We_aq", "We_org", "j_mA_cm2", "FE_ADN_pred"):
        df = df[np.isfinite(df[col]) & (df[col] > 0)].copy()
    df["lwaq"]  = np.log10(df["We_aq"])
    df["lworg"] = np.log10(df["We_org"])

    fe_min, fe_max = df["FE_ADN_pred"].min(), df["FE_ADN_pred"].max()
    pad = 0.05 * (fe_max - fe_min)
    cmap = plt.cm.turbo
    norm = Normalize(vmin=fe_min - pad, vmax=fe_max + pad)
    print(f"\n[{label}]  n={len(df)}  FE_ADN_pred range [{fe_min:.1f}, {fe_max:.1f}] %")

    fig = plt.figure(figsize=(21, 8))
    for k, (gap, glabel) in enumerate(zip(GAP_ORDER, GAP_LABELS)):
        gdf = df[np.isclose(df["gap_mm"], gap, atol=0.01)].copy().reset_index(drop=True)
        ax  = fig.add_subplot(1, 3, k + 1, projection="3d")
        x  = gdf["lwaq"].values
        y  = gdf["lworg"].values
        z  = gdf["j_mA_cm2"].values
        fe = gdf["FE_ADN_pred"].values
        gpr, sc = fit_gpr(x, y, z, fe)

        waq_g  = np.linspace(x.min(), x.max(), N_GRID)
        worg_g = np.linspace(y.min(), y.max(), N_GRID)
        j_g    = np.linspace(z.min(), z.max(), N_GRID)

        WQ, JJ, FE1 = predict_face(gpr, sc, waq_g, j_g, 1, y.max())
        ax.plot_surface(WQ, np.full_like(WQ, y.max()), JJ,
                        facecolors=face_colors(FE1, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)
        WO, JJ, FE2 = predict_face(gpr, sc, worg_g, j_g, 0, x.min())
        ax.plot_surface(np.full_like(WO, x.min()), WO, JJ,
                        facecolors=face_colors(FE2, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)
        WQ, WO, FE3 = predict_face(gpr, sc, waq_g, worg_g, 2, z.max())
        ax.plot_surface(WQ, WO, np.full_like(WQ, z.max()),
                        facecolors=face_colors(FE3, cmap, norm),
                        shade=False, alpha=1.0, linewidth=0)
        ax.scatter(x, y, z, c=fe, cmap=cmap, norm=norm,
                   s=28, edgecolors="k", linewidths=0.4, zorder=5)
        ax.set_title(f"{glabel}\nmodel range: {fe.min():.1f} - {fe.max():.1f} %",
                     fontsize=14, pad=8)
        ax.set_xlabel("We$_\\mathrm{aqueous}$",  fontsize=14, labelpad=14)
        ax.set_ylabel("We$_\\mathrm{organic}$",  fontsize=14, labelpad=14)
        ax.set_zlabel("Current Density (mA cm$^{-2}$)", fontsize=14, labelpad=14)
        xticks = decade_ticks(x); yticks = decade_ticks(y)
        ax.set_xticks(xticks); ax.set_xticklabels([f"$10^{{{t}}}$" for t in xticks], fontsize=12)
        ax.set_yticks(yticks); ax.set_yticklabels([f"$10^{{{t}}}$" for t in yticks], fontsize=12)
        ax.tick_params(axis="z", labelsize=12)
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_zlim(z.min(), z.max())
        ax.set_box_aspect([1, 1, 1])
        ax.view_init(elev=ELEV, azim=AZIM)

    fig.subplots_adjust(left=0.04, right=0.88, wspace=0.20, top=0.86, bottom=0.06)
    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar_ax = fig.add_axes([0.905, 0.18, 0.012, 0.62])
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.ax.tick_params(labelsize=12)
    cb.set_label("FE$_{\\mathrm{ADN}}$  [%]", fontsize=14)
    fig.suptitle(f"Weber-number regime map — {label}\nGPR over modeled "
                 f"FE$_{{\\mathrm{{ADN}}}}$  (n={len(df)} converged rows)",
                 fontsize=14, y=1.00)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_png}")


def main():
    for subdir, label, fname in CASES:
        parity = os.path.join(HERE, "output", subdir, "data", "parity.csv")
        if not os.path.exists(parity):
            print(f"MISSING: {parity}  — run forward_apply_compare.jl first")
            continue
        outdir = os.path.join(HERE, "output", subdir, "data")
        os.makedirs(outdir, exist_ok=True)
        make_one_map(parity, label, os.path.join(outdir, fname))


if __name__ == "__main__":
    main()
