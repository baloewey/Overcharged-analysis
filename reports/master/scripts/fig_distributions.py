#!/usr/bin/env python3
"""
Defect distributions for the master document: time-averaged density maps and wall-distance profiles, same treatment
for circles, Cassini ovals and ellipses, passive and active.

Per run (frames after the first 20% of the run): +1/2 and -1/2 positions from defects.dat (read_with_angle);
  maps     2D histogram, 1-cell bins, per frame and per cell area (defects per cell^2 per frame)
  profile  n_pm(s) = defects per unit area per frame in shells of the distance s to the wall (bin 0.5 cell);
           shell areas from a 0.25-cell grid of points inside the wall (exact for any shape)
  speed    |displacement| / dt of +1/2 matched to the nearest +1/2 of the next frame (<= 3 cells), averaged per shell
Distances use the nearest point of the wall polyline (wall_fourier.wall_polyline / project).
Cached per run in ../data/dist_<name>.npz. Figures in ../figs/.
Usage: fig_distributions.py [--workers 12]
"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUNS = os.path.join(ROOT, "runs")
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle   # noqa: E402
from wall_fourier import wall_polyline, project  # noqa: E402
from scipy.spatial import cKDTree            # noqa: E402

DATA = os.path.join(HERE, "..", "data"); FIGS = os.path.join(HERE, "..", "figs")
LA = 22.21
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
C, K = "cavity/server_main", "cassini"

SETS = {
    "circles": [
        ("passive, R = 1.43 ℓa", f"{C}/main_passive/rho0.06/k06/s01"),
        ("R = 0.57 ℓa", f"{C}/main_act0.07/rho0.15/k06/s01"),
        ("R = 0.86 ℓa", f"{C}/main_act0.07/rho0.10/k06/s01"),
        ("R = 1.43 ℓa", f"{C}/main_act0.07/rho0.06/k06/s01"),
        ("R = 2.15 ℓa", f"{C}/main_act0.07/rho0.04/k06/s01"),
        ("R = 2.87 ℓa", f"{C}/main_act0.07/rho0.03/k06/s01"),
    ],
    "ovals_small": [
        ("passive b/a 2.5", f"{K}/P5_passive_refs/H34_ba2.5_k06"),
        ("passive b/a 2.0", f"{K}/P5_passive_refs/H34_ba2.0_k06"),
        ("passive b/a 1.6", f"{K}/P5_passive_refs/H34_ba1.6_k06"),
        ("passive b/a 1.3", f"{K}/P5_passive_refs/H34_ba1.3_k06"),
        ("b/a 2.5 (aspect 1.17)", f"{K}/P3/H34_ba2.5_k06/s01"),
        ("b/a 2.0 (aspect 1.29)", f"{K}/P3/H34_ba2.0_k06/s01"),
        ("b/a 1.6 (aspect 1.51)", f"{K}/P3/H34_ba1.6_k06/s01"),
        ("b/a 1.3 (aspect 1.94)", f"{K}/P3/H34_ba1.3_k06/s01"),
    ],
    "ovals_large": [
        ("b/a 2.5", f"{K}/W1_size/H48_ba2.5_k06/s01"),
        ("b/a 2.0", f"{K}/W1_size/H48_ba2.0_k06/s01"),
        ("b/a 1.6", f"{K}/W1_size/H48_ba1.6_k06/s01"),
        ("b/a 1.3", f"{K}/W1_size/H48_ba1.3_k06/s01"),
        ("b/a 1.1", f"{K}/P1/H48_ba1.1_k06"),
        ("b/a 1.02 (pinched)", f"{K}/P2/H48_ba1.02_k06"),
        ("passive b/a 1.02", f"{K}/P2/passive_H48_ba1.02_k06"),
        ("passive b/a 1.02, k = 18", f"{K}/P3/passive_H48_ba1.02_k18"),
    ],
    "ellipses": [
        ("passive ε 1.2", f"{K}/P5_passive_refs/ellipse_eps1.2_k06"),
        ("passive ε 1.5", f"{K}/P5_passive_refs/ellipse_eps1.5_k06"),
        ("passive ε 2.0", f"{K}/P5_passive_refs/ellipse_eps2.0_k06"),
        ("ε 1.2", f"{K}/W3_ellipse/eps1.2_k06/s01"),
        ("ε 1.5", f"{K}/W3_ellipse/eps1.5_k06/s01"),
        ("ε 2.0", f"{K}/W3_ellipse/eps2.0_k06/s01"),
    ],
}


def name_of(path):
    return path.replace("/", "__")


def analyse(path):
    out = os.path.join(DATA, f"dist_{name_of(path)}.npz")
    if os.path.exists(out): return out
    run = os.path.join(RUNS, path)
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    W, s, P = wall_polyline(bc)
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); dt = t[1] - t[0]
    ext = np.abs(W).max(0)
    xe = np.arange(-np.ceil(ext[0]) - 1, np.ceil(ext[0]) + 1.01, 1.0); ye = np.arange(-np.ceil(ext[1]) - 1, np.ceil(ext[1]) + 1.01, 1.0)
    path_w = Path(W)
    g = np.mgrid[-ext[0]:ext[0]:0.25, -ext[1]:ext[1]:0.25].reshape(2, -1).T
    g = g[path_w.contains_points(g)]
    dg, _ = project(W, s, g)
    sbins = np.arange(0, dg.max() + 0.5, 0.5)
    shell_area = np.histogram(dg, sbins)[0] * 0.0625
    Hp = np.zeros((len(xe) - 1, len(ye) - 1)); Hm = np.zeros_like(Hp)
    npf = np.zeros(len(sbins) - 1); nmf = np.zeros(len(sbins) - 1); vs = np.zeros(len(sbins) - 1); vn = np.zeros(len(sbins) - 1)
    prev = None
    for f in fr:
        p = f[f[:, 2] > 0, :2]; m = f[f[:, 2] < 0, :2]
        if len(p):
            Hp += np.histogram2d(p[:, 0], p[:, 1], [xe, ye])[0]
            dp, _ = project(W, s, p); npf += np.histogram(dp, sbins)[0]
            if prev is not None and len(prev):
                dd, j = cKDTree(p).query(prev, distance_upper_bound=3.0); ok = np.isfinite(dd)
                if ok.any():
                    dpp, _ = project(W, s, prev[ok])
                    idx = np.clip(np.digitize(dpp, sbins) - 1, 0, len(vs) - 1)
                    np.add.at(vs, idx, dd[ok] / dt); np.add.at(vn, idx, 1)
        if len(m):
            Hm += np.histogram2d(m[:, 0], m[:, 1], [xe, ye])[0]
            dm, _ = project(W, s, m); nmf += np.histogram(dm, sbins)[0]
        prev = p
    nf = len(fr)
    np.savez_compressed(out, Hp=Hp / nf, Hm=Hm / nf, xe=xe, ye=ye, W=W, P=P, k=bc.get("customAnch"),
                        sb=sbins, area=shell_area, np_s=npf / nf / np.maximum(shell_area, 1e-9),
                        nm_s=nmf / nf / np.maximum(shell_area, 1e-9), v_s=vs / np.maximum(vn, 1), vn=vn, nf=nf, dt=dt,
                        Np=npf.sum() / nf, Nm=nmf.sum() / nf)
    return out


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


def gallery(setname, ncol, vmax=None):
    items = SETS[setname]; n = len(items); nrow = int(np.ceil(n / ncol))
    D = [np.load(os.path.join(DATA, f"dist_{name_of(p)}.npz")) for _, p in items]
    vmax = vmax or max(max(d["Hp"].max(), d["Hm"].max()) for d in D)
    hx = max(np.ptp(d["W"][:, 0]) for d in D) / 2 + 2; hy = max(np.ptp(d["W"][:, 1]) for d in D) / 2 + 2
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.9 * nrow * max(hy / hx, 0.45) + 0.5), squeeze=False)
    for ax in axs.flat: ax.axis("off")
    for ax, (lab, _), d in zip(axs.flat, items, D):
        ext = [d["xe"][0], d["xe"][-1], d["ye"][0], d["ye"][-1]]
        for H, cm in ((d["Hp"], "Blues"), (d["Hm"], "Oranges")):
            if H.max() > 0:
                ax.imshow(np.where(H > 0, H, np.nan).T, origin="lower", extent=ext, cmap=cm, interpolation="nearest",
                          norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax))
        ax.plot(np.r_[d["W"][:, 0], d["W"][0, 0]], np.r_[d["W"][:, 1], d["W"][0, 1]], color=INK, lw=0.8)
        ax.set_xlim(-hx, hx); ax.set_ylim(-hy, hy); ax.set_aspect("equal")
        ax.set_title(f"{lab}\n+½: {float(d['Np']):.1f}  −½: {float(d['Nm']):.1f} per frame", fontsize=7.5, color=INK)
    fig.suptitle("Time-averaged defect density (log scale, defects per cell² per frame): +½ blue, −½ orange; "
                 f"same scale in every panel (max {vmax:.2g}); k = 6 unless noted", fontsize=8, color=INK2)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, f"maps_{setname}.png"), dpi=170); plt.close(fig)


def profiles(setname, title):
    items = SETS[setname]
    D = [np.load(os.path.join(DATA, f"dist_{name_of(p)}.npz")) for _, p in items]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(items)))
    for (lab, _), d, c in zip(items, D, cols):
        sc = 0.5 * (d["sb"][1:] + d["sb"][:-1]) / LA
        ls = "--" if lab.startswith("passive") else "-"
        ok = d["area"] > 20
        axs[0].plot(sc[ok], d["np_s"][ok], ls, color=c, lw=1.5, label=lab)
        axs[1].plot(sc[ok], d["nm_s"][ok], ls, color=c, lw=1.5)
        okv = ok & (d["vn"] > 200)
        axs[2].plot(sc[okv], d["v_s"][okv], ls, color=c, lw=1.5)
    for ax, yl in zip(axs, ("+½ density (per cell² per frame)", "−½ density (per cell² per frame)", "+½ speed (cells per time unit)")):
        ax.set_xlabel("distance from the wall s / ℓa", fontsize=9, color=INK); ax.set_ylabel(yl, fontsize=9, color=INK)
        ax.set_yscale("log" if "density" in yl else "linear"); style(ax); ax.set_xlim(0, None)
    axs[0].legend(frameon=False, fontsize=6.5)
    fig.suptitle(title + " (dashed: passive). Speed from nearest-neighbour matching of +½ between frames 2 time units apart; "
                 "passive speed = detection jitter", fontsize=8, color=INK2)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, f"profiles_{setname}.png"), dpi=170); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
    paths = [p for v in SETS.values() for _, p in v]
    todo = [p for p in paths if os.path.exists(os.path.join(RUNS, p, "DONE"))]
    skipped = sorted(set(paths) - set(todo))
    for p in skipped: print("not finished, skipped:", p)
    with Pool(a.workers) as pool:
        for o in pool.imap_unordered(analyse, todo): print("ok", os.path.basename(o), flush=True)
    for setname, ncol, title in (("circles", 3, "Circles, k = 6"), ("ovals_small", 4, "Cassini ovals, half-length 1.53 ℓa, k = 6"),
                                 ("ovals_large", 4, "Cassini ovals, half-length 2.16 ℓa, k = 6"), ("ellipses", 3, "Ellipses, area 2500, k = 6")):
        if all(p in todo for _, p in SETS[setname]):
            gallery(setname, ncol); profiles(setname, title); print("figures", setname)
        else:
            print("set incomplete, no figure yet:", setname)


if __name__ == "__main__":
    main()
