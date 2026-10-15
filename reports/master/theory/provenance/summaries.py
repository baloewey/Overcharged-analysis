#!/usr/bin/env python3
"""Summary tables for numbers quoted in the papers that had no output file (PRL number audit, user 2026-10-15).
  swirl_summary.txt     passive MPCD cavities (k = 6) under an imposed vortex (runs/swirl_threshold, runs/swirl_tumble): wall +1/2 count (within 4 cells),
                        bulk net charge, tilt coherence C = |<exp(i tau)>| of the wall +1/2 (tau = tail angle - outward normal); second half of each run.
  neglayer_summary.txt  negative winding k = -6 (runs/neglayer, area 2500, seeded starts): wall -1/2 count, bulk +1/2 and -1/2, t >= 2500 t.u.
  perlength_qbins.txt   wall law from energetics/perlength.npy (shape sweep): per q bin, aspect >= 1.5: bulk charge per unit wall length c and released fraction Q_b/k.
Usage: summaries.py"""
import os, sys, glob, json, numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa


def per_run(run, tmin=None, frac=0.5):
    bc = json.load(open(run + "/input.json"))["BC"][0]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
    sel = (t >= tmin) if tmin is not None else (t >= t[0] + frac * (t[-1] - t[0]))
    wp, wm, bp, bm, co = [], [], [], [], []
    for f in [fr[i] for i in np.flatnonzero(sel)]:
        if not len(f): continue
        d, j = tree.query(f[:, :2]); w = d < 4
        wp.append(((f[:, 2] > 0) & w).sum()); wm.append(((f[:, 2] < 0) & w).sum()); bp.append(((f[:, 2] > 0) & ~w).sum()); bm.append(((f[:, 2] < 0) & ~w).sum())
        ww = w & (f[:, 2] > 0)
        if ww.any(): co.append(abs(np.mean(np.exp(1j * (f[ww, 3] - np.angle(Nn[j[ww]]))))))
    return bc.get("customAnch"), np.mean(wp), np.mean(wm), 0.5 * (np.mean(bp) - np.mean(bm)), np.mean(bp), np.mean(bm), np.mean(co) if co else np.nan, len(wp)


if __name__ == "__main__":
    L = ["IMPOSED VORTEX on passive MPCD cavities (k = 6; second half of each run; wall = within 4 cells)",
         "run | k | wall +1/2 | wall -1/2 | bulk net charge | tilt coherence C (wall +1/2) | frames"]
    for run in sorted(glob.glob(ROOT + "/runs/swirl_threshold/passive_*") + glob.glob(ROOT + "/runs/swirl_tumble/passive_*")):
        k, wp, wm, qb, bp, bm, co, n = per_run(run); L.append(f"{os.path.relpath(run, ROOT + '/runs'):40s} | {k} | {wp:5.2f} | {wm:4.2f} | {qb:+.2f} | {co:.3f} | {n}")
    open(os.path.join(H, "swirl_summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    L = ["NEGATIVE WINDING k = -6 (runs/neglayer; area 2500; t >= 2500 t.u.; wall = within 4 cells). Continuum prediction 2|k|+2 = 14 wall -1/2.",
         "run | wall -1/2 | wall +1/2 | bulk +1/2 | bulk -1/2 | bulk net | frames"]
    for run in sorted(glob.glob(ROOT + "/runs/neglayer/*")):
        if not os.path.exists(run + "/input.json"): continue
        k, wp, wm, qb, bp, bm, co, n = per_run(run, tmin=2500); L.append(f"{os.path.basename(run):28s} | {wm:5.2f} | {wp:4.2f} | {bp:4.2f} | {bm:4.2f} | {qb:+.2f} | {n}")
    open(os.path.join(H, "neglayer_summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    A = np.load(ROOT + "/reports/master/theory/energetics/perlength.npy", allow_pickle=True)
    fam, asp, req, P, k, nw, qb = (np.array([a[i] for a in A]) for i in range(7)); asp, P, k, qb = (x.astype(float) for x in (asp, P, k, qb)); q = k / P
    L = ["WALL LAW per winding-density bin (energetics/perlength.npy, shape sweep, 3 seeds; aspect >= 1.5): c = Q_b/P (bulk charge per cell of wall length), released fraction Q_b/k",
         "q bin | runs | q mean | c mean +- s.e. | Q_b/k mean +- s.e."]
    for lo, hi in ((0.015, 0.025), (0.025, 0.037), (0.037, 0.052), (0.052, 0.065), (0.065, 0.08)):
        m = (q >= lo) & (q < hi) & (asp >= 1.5)
        if not m.any(): continue
        c = qb[m] / P[m]; f = qb[m] / k[m]
        L.append(f"[{lo:.3f}, {hi:.3f}) | {m.sum():3d} | {q[m].mean():.4f} | {c.mean():.4f} +- {c.std(ddof=1)/np.sqrt(m.sum()):.4f} | {f.mean():.3f} +- {f.std(ddof=1)/np.sqrt(m.sum()):.3f}")
    open(os.path.join(H, "perlength_qbins.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
