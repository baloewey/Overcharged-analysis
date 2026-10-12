#!/usr/bin/env python3
"""
Passive circles (activity 0) across the main sweep: where the imposed defects sit (master document).

Runs: runs/cavity/server_main/main_passive (+k) and main_passive_neg (-k), 1 seed each, 10k warm-up + 30k steps.
Per run (frames after the first 20%): all defects; s = R - r (distance from the wall).
  N_+, N_- per frame; mean and 95th percentile of s for the imposed-sign defects; q_bulk(s0) = net charge of defects
  with s > s0 per frame (s0 = 5 cells); frames with exactly 2|k| imposed-sign defects and no other defect.
Figure ../figs/passive_circles.png: (a) maps for k = 1, 2, 3, 6, 12 at rho = 0.06 and k = -6; (b) distribution of s for
k = 1..12 at rho = 0.06; (c) 95th percentile of s over (rho, R/l_a); (d) mean distance of the imposed defects from the wall vs the imposed spacing 1/rho.
Table ../data/passive_circles.tsv.
"""
import argparse, glob, json, os, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "runs"))
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle   # noqa: E402
DATA = os.path.join(HERE, "..", "data"); FIGS = os.path.join(HERE, "..", "figs")
LA = 22.21; S0 = 5.0
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"


def one(run):
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    k = int(bc["customAnch"]); R = bc["R"]; sg = 1 if k > 0 else -1
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2)
    sA, Np, Nm, qb, clean = [], [], [], [], 0
    pos = []
    for f in fr:
        r = np.hypot(f[:, 0], f[:, 1]); s = R - r; q = np.sign(f[:, 2])
        Np.append((q > 0).sum()); Nm.append((q < 0).sum())
        sA.append(s[q == sg]); qb.append(0.5 * q[s > S0].sum())
        clean += ((q == sg).sum() == 2 * abs(k)) and ((q == -sg).sum() == 0)
        pos.append(f[:, :3])
    sA = np.concatenate(sA) if sA else np.array([])
    return dict(run=run, k=k, R=R, rho=abs(k) / (np.pi * R), Np=np.mean(Np), Nm=np.mean(Nm), s_mean=sA.mean(), s95=np.percentile(sA, 95),
                qb=np.mean(qb), clean=clean / len(fr), s_all=sA, pos=np.concatenate(pos))


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    base = os.path.join(RUNS, "cavity", "server_main")
    runs = sorted(glob.glob(os.path.join(base, "main_passive*", "rho*", "k*", "s01")))
    runs = [r for r in runs if os.path.exists(os.path.join(r, "DONE"))]
    with Pool(a.workers) as pool: res = pool.map(one, runs)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "passive_circles.tsv"), "w") as f:
        f.write("set\trho\tk\tR\tR/la\tN+\tN-\ts_mean\ts95\tq_bulk(s>5)\tclean_frames\n")
        for r in res:
            st = os.path.relpath(r["run"], base).split("/")[0]
            f.write(f"{st}\t{r['rho']:.3f}\t{r['k']}\t{r['R']:.2f}\t{r['R'] / LA:.2f}\t{r['Np']:.2f}\t{r['Nm']:.2f}\t{r['s_mean']:.2f}\t{r['s95']:.2f}\t{r['qb']:.3f}\t{r['clean']:.3f}\n")
    pos = [r for r in res if r["k"] > 0]; neg = [r for r in res if r["k"] < 0]
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.15])
    pick = [(1, 0.06), (2, 0.06), (3, 0.06), (6, 0.06), (12, 0.06), (-6, 0.06)]
    for i, (k, rho) in enumerate(pick):
        ax = fig.add_subplot(gs[0, i])
        r = min((x for x in res if x["k"] == k), key=lambda x: abs(x["rho"] - rho))
        p = r["pos"]; R = r["R"]; e = np.arange(-R - 1, R + 1.01, 0.5)
        for sgn, cm in ((1, "Blues"), (-1, "Oranges")):
            q = p[np.sign(p[:, 2]) == sgn]
            if len(q):
                H = np.histogram2d(q[:, 0], q[:, 1], [e, e])[0]
                ax.imshow(np.where(H > 0, H, np.nan).T, origin="lower", extent=[e[0], e[-1], e[0], e[-1]], cmap=cm, norm=LogNorm(1, H.max()), interpolation="nearest")
        th = np.linspace(0, 2 * np.pi, 300); ax.plot(R * np.cos(th), R * np.sin(th), color=INK, lw=0.8)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"k = {k}, R = {R:.1f}\n({R / LA:.2f} ℓa, ρ = {r['rho']:.2f})", fontsize=8, color=INK)
    # (b) s distributions at rho = 0.06
    ax = fig.add_subplot(gs[1, 0:2])
    sel = sorted([x for x in pos if abs(x["rho"] - 0.06) < 0.003], key=lambda x: x["k"])
    for x in sel:
        h, e = np.histogram(x["s_all"], np.arange(0, 20.5, 0.5), density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color=plt.cm.viridis(x["k"] / 13), lw=1.4, label=f"k = {x['k']}")
    ax.set_xlabel("distance from the wall s (cells)", fontsize=9, color=INK); ax.set_ylabel("probability density", fontsize=9, color=INK)
    ax.set_title("(b) imposed defects, ρ = 0.06", loc="left", fontsize=10, color=INK); ax.legend(frameon=False, fontsize=6.5, ncol=2); style(ax)
    # (c) s95 map
    ax = fig.add_subplot(gs[1, 2:4])
    for xs, mk, off in ((pos, "o", 0), (neg, "s", 0.004)):
        sc = ax.scatter([x["rho"] + off for x in xs], [x["R"] / LA for x in xs], c=[x["s95"] for x in xs], cmap="viridis", vmin=1, vmax=6,
                        marker=mk, s=28, edgecolor="none")
    fig.colorbar(sc, ax=ax, fraction=0.05)
    ax.set_xlabel("ρ = 2|k|/P", fontsize=9, color=INK); ax.set_ylabel("R / ℓa", fontsize=9, color=INK)
    ax.set_title("(c) 95th percentile of s (cells); ○ +k, □ −k", loc="left", fontsize=10, color=INK); style(ax)
    # (d) passive layer depth vs imposed spacing
    ax = fig.add_subplot(gs[1, 4:6])
    for xs, mk, lab in ((pos, "o", "+k"), (neg, "s", "−k")):
        ax.scatter([1 / x["rho"] for x in xs], [x["s_mean"] for x in xs], c=[x["R"] / LA for x in xs], cmap="magma_r", vmin=0, vmax=3,
                   marker=mk, s=24, edgecolor="none", label=lab)
    sp = np.array([1 / x["rho"] for x in res if x["rho"] <= 0.151]); sm = np.array([x["s_mean"] for x in res if x["rho"] <= 0.151])
    c1, c0 = np.polyfit(sp, sm, 1); xx = np.linspace(0, 35, 10)
    ax.plot(xx, c0 + c1 * xx, color=INK2, lw=1, ls="--", label=f"fit ρ ≤ 0.15: s = {c0:.2f} + {c1:.3f}/ρ")
    ax.set_xlabel("imposed spacing 1/ρ = P/(2|k|) (cells)", fontsize=9, color=INK); ax.set_ylabel("mean distance from the wall (cells)", fontsize=9, color=INK)
    ax.set_title("(d) passive layer depth vs spacing (colour: R/ℓa)", loc="left", fontsize=10, color=INK); ax.legend(frameon=False, fontsize=7); style(ax)
    print(f"depth fit: s = {c0:.3f} + {c1:.4f}/rho")
    fig.suptitle("Passive circles (activity 0), main sweep, 1 seed: defect maps (log density; +½ blue, −½ orange), distance from the wall, bulk charge",
                 fontsize=9, color=INK2)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "passive_circles.png"), dpi=170)
    q = np.array([x["qb"] for x in res]); cl = np.array([x["clean"] for x in res])
    print(f"{len(res)} runs; |q_bulk| max {np.abs(q).max():.3f} (k: {[x['k'] for x in res if abs(x['qb']) > 0.05]}); clean frames median {np.median(cl):.3f}, min {cl.min():.3f}")
    for x in sorted(sel, key=lambda x: x["k"]): print(f"k {x['k']:2d} R {x['R']:.1f} s_mean {x['s_mean']:.2f} s95 {x['s95']:.2f} N+ {x['Np']:.2f} N- {x['Nm']:.2f} qb {x['qb']:.3f}")


if __name__ == "__main__":
    main()
