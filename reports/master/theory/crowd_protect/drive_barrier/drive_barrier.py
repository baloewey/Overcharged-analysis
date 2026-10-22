#!/usr/bin/env python3
"""Drive x barrier test (user, 2026-10-22). Hypothesis: escape = drive (geometry: coherence loss turns the flow against the layer) x barrier (elastic pull
of the layer, set by the wall density). Shape sweep (675 runs), +1/2 tracked with NEW v2 (tracker_validation/trackers.track_new, gate 9).
Starts: wall residents (d < 4 for the previous 20 t.u.) every 10 t.u.; for each, the deepest point reached in the next 20 t.u., x_max.
  S(x) = P(x_max > x), x = 1.5..10 cells (the crowded-layer force peak sits 1.5-2 cells beyond the resident depth, inside the 4-cell band).  Attempt A = S(4) (leaves the band).  Completion c6 = S(6)/S(4), c8 = S(8)/S(4).
  Prediction 1: A follows coherence (drive), not rho; c6, c8 follow rho (barrier), not coherence.
  Prediction 2 (effective temperature): ln S(x) = a - U(x)/T_a with U(x) the one-constant elastic cost of dragging a wall +1/2 from its resident depth d0
  to x (force_vs_N.npz, circle P 200, interpolated to N = rho x 200); T_a depends on the drive (coherence), not on rho, and -ln S/U collapses across rho.
-> drive_barrier.csv (per run), drive_barrier.txt (summary)"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__)); XG = np.arange(1.5, 10.01, 0.5)
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr)
        T = [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]; xm = []; d0s = []
        for tr in T:
            nn = np.array(tr["n"]); d = dist(np.array(tr["xy"])); pos = dict(zip(nn, d))
            for n in nn:
                if n % 5 or n + 10 >= nf or n < 10: continue
                hist = [pos.get(m) for m in range(n - 10, n + 1)]
                if any(h is None or h >= 4 for h in hist): continue
                fut = [pos[m] for m in range(n + 1, n + 11) if m in pos]
                if not fut: continue
                xm.append(max(fut)); d0s.append(pos[n])
        xm = np.array(xm); S = [(xm > x).mean() if len(xm) else np.nan for x in XG]
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), n=len(xm), d0=float(np.median(d0s)) if d0s else np.nan, **{f"S{x:g}": s for x, s in zip(XG, S)})
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "n", "d0"] + [f"S{x:g}" for x in XG] + ["err"]
    with open(os.path.join(H, "drive_barrier.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
