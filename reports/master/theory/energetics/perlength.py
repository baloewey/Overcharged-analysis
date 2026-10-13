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
    fig, ax = plt.subplots(1, 3, figsize=(14, 4)); cm = {1.0: "#1f5fa8", 1.5: "#7a3b8f", 2.0: "#c2452d"}
    for fam, mk in (("circle", "*"), ("ellipse", "o"), ("oval", "s")):
        for r in (1.0, 1.5, 2.0):
            S = [a for a in A if a[0] == fam and a[2] == r]
            if not S: continue
            asp = np.array([a[1] for a in S]); ax[0].scatter(asp + 0.01 * (hash(fam) % 3), [a[6] / a[3] for a in S], s=12, marker=mk, color=cm[r], alpha=0.6, label=f"{fam} R_eq {r}" if fam != "circle" or r == 1.0 else None)
            ax[1].scatter(asp, [a[6] / a[4] for a in S], s=12, marker=mk, color=cm[r], alpha=0.6); ax[2].scatter([a[4] / a[3] for a in S], [a[6] / a[3] for a in S], s=12, marker=mk, color=cm[r], alpha=0.6)
    ax[0].axhspan(0.012, 0.0145, color="0.85", label="channel (q = 0.031)"); ax[2].axhspan(0.012, 0.0145, color="0.85"); ax[2].axvline(0.031, color="0.6", ls=":")
    ax[0].set_xlabel("aspect ratio"); ax[0].set_ylabel("bulk charge per unit wall length"); ax[1].set_xlabel("aspect ratio"); ax[1].set_ylabel("released fraction Q_b / k")
    ax[2].set_xlabel("winding density q = k/P"); ax[2].set_ylabel("bulk charge per unit wall length"); ax[0].legend(fontsize=6, frameon=False, ncol=2)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_perlength.png"), dpi=200)
