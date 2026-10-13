#!/usr/bin/env python3
"""Charge per unit wall length vs aspect ratio (user, 2026-10-13). Shape sweep (ellipses, ovals, circles), second half, every 5th frame.
Per run: perimeter P, winding k, winding density q = k/P, wall +1/2 (< 4 cells), bulk net charge Q_b (>= 4 cells); Q_b/P, wall +1/2 per length, released fraction Q_b/k.
Channel reference (runs/channel, ellipse winding density): Q_b per wall length 0.012-0.014, released fraction 0.39-0.46. -> perlength.txt, fig_perlength.png"""
import os, sys, csv, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
        t, fr = En.read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); nw, qb = [], []
        for f in fr[::5]:
            if not len(f): continue
            d = tree.query(f[:, :2])[0]; nw.append(((f[:, 2] > 0) & (d < 4)).sum()); qb.append(0.5 * np.sign(f[d >= 4, 2]).sum())
        return (rel, P, k, np.mean(nw), np.mean(qb))
    except Exception as ex:
        return (rel, None, str(ex))
if __name__ == "__main__":
    if os.environ.get("REPLOT"):
        A = [tuple(a) for a in np.load(os.path.join(H, "perlength.npy"), allow_pickle=True)]
    else:
        man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
        rels = [r for r in man if os.path.exists(os.path.join(SW, r, "DONE"))]
        with Pool(8) as pool: res = [r for r in pool.map(one, rels, chunksize=4) if r[1] is not None]
        A = [(man[r[0]]["family"], float(man[r[0]]["aspect"]), float(man[r[0]]["Req_la"]), r[1], r[2], r[3], r[4]) for r in res]
        np.save(os.path.join(H, "perlength.npy"), np.array(A, dtype=object), allow_pickle=True)
    L = ["family aspect | runs | winding density q = k/P (range) | wall +1/2 per length | bulk charge per length | released fraction Q_b/k | (by R_eq 1.0 / 1.5 / 2.0: Q_b per length)"]
    for fam in ("circle", "ellipse", "oval"):
        for asp in sorted(set(round(a[1], 2) for a in A if a[0] == fam)):
            S = [a for a in A if a[0] == fam and round(a[1], 2) == asp]
            q = np.array([a[4] / a[3] for a in S]); wl = np.array([a[5] / a[3] for a in S]); bl = np.array([a[6] / a[3] for a in S]); fr_ = np.array([a[6] / a[4] for a in S])
            byR = " / ".join(f"{np.mean([a[6] / a[3] for a in S if a[2] == r]):.4f}" for r in (1.0, 1.5, 2.0) if any(a[2] == r for a in S))
            L.append(f"{fam:7s} {asp:4.2f} | {len(S)} | {q.min():.3f}-{q.max():.3f} | {wl.mean():.4f} | {bl.mean():.4f} ± {bl.std() / np.sqrt(len(bl)):.4f} | {fr_.mean():.2f} | {byR}")
    print("\n".join(L)); open(os.path.join(H, "perlength.txt"), "w").write("\n".join(L) + "\n")
    # Colour = winding density (the variable that sets the plateau); marker = family. Channel (runs/channel, q 0.031, W 20/40/60) drawn as the
    # flat limit: an extra x position "channel" in (a, b) and a point at q 0.031 in (c), in the colour of its q bin. (Earlier version: grey bands, unclear.)
    CH_Q, CH_QB, CH_QBE, CH_F, CH_FE = 0.031, 0.0132, 0.0004, 0.42, 0.01
    qbins = [(0.0, 0.025, "q ≈ 0.02"), (0.025, 0.037, "q ≈ 0.03"), (0.037, 0.052, "q ≈ 0.045"), (0.052, 0.065, "q ≈ 0.06"), (0.065, 1, "q ≈ 0.07")]
    qcol = [matplotlib.colormaps["viridis"](v) for v in (0.0, 0.25, 0.5, 0.72, 0.9)]  # ordered: low q dark, high q light
    def qi(q): return next(i for i, (lo, hi, _) in enumerate(qbins) if lo <= q < hi)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2)); off = {"circle": 0.0, "ellipse": -0.012, "oval": 0.012}; XC = 2.85
    for fam, mk in (("circle", "*"), ("ellipse", "o"), ("oval", "s")):
        S = [a for a in A if a[0] == fam]
        q = np.array([a[4] / a[3] for a in S]); c = [qcol[qi(v)] for v in q]; asp = np.array([a[1] for a in S]) + off[fam]
        ax[0].scatter(asp, [a[6] / a[3] for a in S], s=14 if mk != "*" else 30, marker=mk, c=c, alpha=0.7, lw=0)
        ax[1].scatter(asp, [a[6] / a[4] for a in S], s=14 if mk != "*" else 30, marker=mk, c=c, alpha=0.7, lw=0)
        ax[2].scatter(q, [a[6] / a[3] for a in S], s=14 if mk != "*" else 30, marker=mk, c=["0.55" if a[1] < 1.3 else "k" for a in S], alpha=0.6, lw=0)
    cc = qcol[qi(CH_Q)]
    ax[0].errorbar([XC], [CH_QB], yerr=[CH_QBE], fmt="D", ms=8, color=cc, mec="k", capsize=3)
    ax[1].errorbar([XC], [CH_F], yerr=[CH_FE], fmt="D", ms=8, color=cc, mec="k", capsize=3)
    ax[2].errorbar([CH_Q], [CH_QB], yerr=[CH_QBE], fmt="D", ms=8, color="w", mec="k", ecolor="k", capsize=3, zorder=5)
    for a_ in ax[:2]:
        a_.axvline(2.68, color="0.7", lw=0.8, ls=":"); a_.set_xticks([1.0, 1.4, 1.8, 2.2, 2.6 - 0.1, XC]); a_.set_xticklabels(["1.0", "1.4", "1.8", "2.2", "2.5", "channel"]); a_.set_xlim(0.92, 3.0); a_.set_xlabel("aspect ratio")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=qcol[i], label=lab) for i, (_, _, lab) in enumerate(qbins)]
    h += [Line2D([], [], marker=m, ls="", color="0.4", label=l) for m, l in (("*", "circle"), ("o", "ellipse"), ("s", "oval"))]
    h += [Line2D([], [], marker="D", ls="", color="w", mec="k", label="channel, q = 0.031")]
    fig.legend(handles=h, fontsize=8, frameon=False, ncol=9, loc="lower center", bbox_to_anchor=(0.5, 0.0))
    ax[0].set_ylabel("bulk charge per unit wall length"); ax[1].set_ylabel("released fraction $Q_b/k$")
    ax[2].set_xlabel("winding density q = k/P"); ax[2].set_ylabel("bulk charge per unit wall length")
    ax[2].legend(handles=[Line2D([], [], marker="o", ls="", color="k", label="aspect ≥ 1.3"), Line2D([], [], marker="o", ls="", color="0.55", label="aspect < 1.3"),
                          Line2D([], [], marker="D", ls="", color="w", mec="k", label="channel")], fontsize=7, frameon=False)
    for i, t in enumerate(("(a) rises with elongation, then saturates", "(b) fraction of the imposed charge in the bulk", "(c) the plateau is set by the winding density")):
        ax[i].set_title(t, fontsize=9)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(os.path.join(H, "fig_perlength.png"), dpi=200)
