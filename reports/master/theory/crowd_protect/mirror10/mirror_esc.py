#!/usr/bin/env python3
"""Mirror control (-k) at 10 seeds with the validated tracker (27 Oct list, item 3). runs/shapes_mirror (2250 runs: the 675-run shape grid with k -> -k, 10 seeds).
For -k the wall layer holds -1/2 (2|k| + 2 expected). Per run (second half): wall -1/2 count (d < 4) and its deviation from 2|k|+2; escape probability of
wall -1/2 (shape_runs definition: start d < 4 every 10 t.u., deeper than 6 cells 20 t.u. later), tracks from tracker_validation/trackers.track_new (gate 9),
-1/2 sign. Wall -1/2 sit 2-4 cells deep, so the 6-cell criterion is reached by jitter: also p_strict = resident (d < 4 for the
previous 20 t.u.) and deeper than 10 cells at +20 t.u. -> mirror_esc.csv"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr); k = abs(int(json.load(open(run + "/input.json"))["BC"][0]["customAnch"]))
        nw = np.mean([((dist(f[:, :2]) < 4) & (f[:, 2] < 0)).sum() for f in fr if len(f)])
        T = [tr for tr in track_new(fr, 9.0) if tr["sg"] < 0]; st = es = st2 = es2 = 0
        for tr in T:
            pos = dict(zip(tr["n"], dist(np.array(tr["xy"]))))
            for n, dv in pos.items():
                if n % 5 or n + 10 >= nf or dv >= 4: continue
                st += 1; e = pos.get(n + 10); es += e is not None and e > 6
                if all(pos.get(m) is not None and pos[m] < 4 for m in range(n - 10, n)):      # strict: resident 20 t.u., deeper than 10 cells
                    st2 += 1; es2 += e is not None and e > 10
        return dict(run=os.path.relpath(run, ROOT + "/shapes_mirror"), k=k, nwall=nw, dev=nw - (2 * k + 2), p_esc=es / max(st, 1), starts=st, p_strict=es2 / max(st2, 1))
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_mirror"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_mirror/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(one, runs, chunksize=4)
    with open("mirror_esc.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["run", "k", "nwall", "dev", "p_esc", "starts", "p_strict", "err"], extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
