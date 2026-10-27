#!/usr/bin/env python3
"""+k counterpart of mirror_esc.py's strict escape (resident d < 4 for 20 t.u., deeper than 10 cells 20 t.u. later), +1/2, runs/shapes_server (675 runs, 3 seeds),
for the +k / -k escape ratio of the mirror control. -> plus_strict.csv"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr); st = es = 0
        for tr in [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]:
            pos = dict(zip(tr["n"], dist(np.array(tr["xy"]))))
            for n, dv in pos.items():
                if n % 5 or n + 10 >= nf or dv >= 4 or not all(pos.get(m) is not None and pos[m] < 4 for m in range(n - 10, n)): continue
                st += 1; e = pos.get(n + 10); es += e is not None and e > 10
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), p_strict=es / max(st, 1), starts=st)
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 10))) as p: res = p.map(one, runs, chunksize=2)
    with open("plus_strict.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["run", "p_strict", "starts", "err"], extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
