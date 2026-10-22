#!/usr/bin/env python3
"""Recheck of the shape_runs.csv escape probability with the predictive tracker (user, 2026-10-22). Definition as shape_runs: every wall +1/2 (d < 4)
at start frames every 10 t.u. is followed 20 t.u.; escape = found deeper than 6 cells at +20 t.u. Probability per wall defect per 20 t.u. Here the
follow-up uses exchange_count2.tracks (all frames, 2 t.u.; prediction + Hungarian, gate 9, 2-frame gaps; a track lost before +20 t.u. counts as no escape,
as in the original). Also the lagged outcome of the per-event analysis (in the wall at +40 t.u., deeper than 6 at +60). -> shape_esc2.csv"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
import exchange_count as X, exchange_count2 as X2
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr)
        T = [tr for tr in X2.tracks(fr, 9.0) if tr["sg"] > 0]; st = es = la = 0
        for tr in T:
            nn = np.array(tr["n"]); d = dist(np.array(tr["xy"])); pos = {n: dv for n, dv in zip(nn, d)}
            for n, dv in zip(nn, d):
                if n % 5 or dv >= 4 or n + 10 >= nf: continue
                st += 1; d10 = pos.get(n + 10)
                if d10 is not None and d10 > 6: es += 1
                d20, d30 = pos.get(n + 20), pos.get(n + 30)
                if d20 is not None and d30 is not None and d20 < 4 and d30 > 6: la += 1
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), starts=st, esc=es, late=la, p_esc=es / max(st, 1), p_late=la / max(st, 1))
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 6))) as p: res = p.map(one, runs, chunksize=2)
    with open(os.path.join(H, "shape_esc2.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["run", "starts", "esc", "late", "p_esc", "p_late", "err"], extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
