#!/usr/bin/env python3
"""Item 6 (27 Oct): why do circles respond weakly to crowding? Circle runs of the shape sweep (135), NEW v2 tracks. Excursion profile S(x) (deepest point within
20 t.u. of resident wall +1/2, as drive_barrier.py), split by the layer's occupancy at the start frame: FULL (wall +1/2 count = 2k-2), OVER-FULL (> 2k-2),
DEPLETED (< 2k-2). Near-layer active temperature per class from ln S = a - U/T over x in [d0+0.5, d0+2.5] (relaxed-neighbour U), pooled per rho level.
Expectation if the circle exception is composition: FULL rings follow the elongated-cavity law (T independent of rho); over-full shedding differs. -> circle_occupancy.txt"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT; XG = np.arange(1.5, 10.01, 0.5)
def one(run):
    kind, Q, dist, Lw, info = X.geometry(run); k = int(info["k"]); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr)
    cnt = [((dist(f[:, :2]) < 4) & (f[:, 2] > 0)).sum() if len(f) else 0 for f in fr]; xm = {"full": [], "over": [], "dep": []}; d0 = []
    for tr in [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]:
        pos = dict(zip(tr["n"], dist(np.array(tr["xy"]))))
        for n in tr["n"]:
            if n % 5 or n + 10 >= nf or n < 10 or any(pos.get(m) is None or pos[m] >= 4 for m in range(n - 10, n + 1)): continue
            fut = [pos[m] for m in range(n + 1, n + 11) if m in pos]
            if not fut: continue
            c = "full" if cnt[n] == 2 * k - 2 else ("over" if cnt[n] > 2 * k - 2 else "dep"); xm[c].append(max(fut)); d0.append(pos[n])
    return os.path.relpath(run, ROOT + "/shapes_server"), info["rho"], np.median(d0) if d0 else np.nan, {c: np.array(v) for c, v in xm.items()}
if __name__ == "__main__":
    Z = np.load("../drag_relaxed.npz"); NS = np.array([8, 12, 18, 24, 28]); DEL = Z["DEL"]; EN = np.array([Z[f"relaxed_N{N}"] for N in NS])
    U = lambda x, d0, rho: np.interp(np.clip(x - d0, 0, DEL[-1]), DEL, np.array([np.interp(np.clip(rho * 200, 8, 28), NS, EN[:, j]) for j in range(len(DEL))]))
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/circle/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 10))) as p: res = p.map(one, runs, chunksize=1)
    LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); L = ["CIRCLES BY OCCUPANCY (item 6): per rho level, pooled over runs: starts share | S(d0+1) | S(d0+3) | near-layer T (kBT) | starts"]
    for i, r0 in enumerate(LV):
        g = [r for r in res if np.argmin(abs(r[1] - LV)) == i]; d0 = np.nanmedian([r[2] for r in g]); tot = sum(sum(len(v) for v in r[3].values()) for r in g); row = []
        for c in ("full", "over", "dep"):
            x = np.concatenate([r[3][c] for r in g]) if g else np.zeros(0)
            if len(x) < 50: row.append(f"{c}: n {len(x)}"); continue
            S = np.array([(x > xx).mean() for xx in XG]); ok = (XG >= d0 + 0.5) & (XG <= d0 + 2.5) & (S * len(x) >= 5) & (S > 0); T = np.nan
            if ok.sum() >= 3:
                b = np.polyfit([U(xx, d0, r0) for xx in XG[ok]], np.log(S[ok]), 1); T = -1 / b[0] if b[0] < 0 else np.nan
            row.append(f"{c}: {len(x) / tot:.2f} | {np.interp(d0 + 1, XG, S):.4f} | {np.interp(d0 + 3, XG, S):.4f} | {T:5.1f} | {len(x)}")
        L.append(f"rho {r0:.2f}:  " + "  ||  ".join(row))
    open("circle_occupancy.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
