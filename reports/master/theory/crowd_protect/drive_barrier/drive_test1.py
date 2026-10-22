#!/usr/bin/env python3
"""Test 1 of the drive x barrier hypothesis (user, 2026-10-22: 'run test 1 tonight'). Question: is the extra suppression of near-layer excursions in crowded
layers explained by a weaker drive? Per run (shape sweep, NEW v2 tracks): resident wall +1/2 (d < 4 for the previous 20 t.u.); their depth MSD at lags 1 and 2
frames (2 and 4 t.u.); drive = depth diffusivity D_r = (MSD(2) - MSD(1)) / (2 * 2 t.u.) (removes the static half-cell position noise; overdamped picture:
an active temperature T_a = gamma D_r). Also <v_r^2> at lag 1 for reference.
Barrier side: near-layer active temperature T_near from the excursion profile S(x) of drive_barrier.csv, fitted only from x = d0 + 0.5 to d0 + 2.5
(the deep tail is kinetic), against the relaxed-neighbour drag energy U(x). PASS if T_near / D_r is the same across rho and coherence classes.
-> drive_test1.csv, drive_test1.txt"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); dt = t[1] - t[0]
        T = [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]; m1 = []; m2 = []
        for tr in T:
            nn = np.array(tr["n"]); d = dist(np.array(tr["xy"])); pos = dict(zip(nn, d))
            for n in nn:
                w = [pos.get(m) for m in range(n - 10, n + 3)]
                if any(v is None for v in w) or max(w[:11]) >= 4: continue
                m1.append((pos[n + 1] - pos[n]) ** 2); m2.append((pos[n + 2] - pos[n]) ** 2)
        M1, M2 = np.mean(m1), np.mean(m2)
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), n=len(m1), msd1=M1, msd2=M2, D_r=(M2 - M1) / (2 * dt), v2=M1 / dt ** 2)
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    with open(os.path.join(H, "drive_test1.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["run", "n", "msd1", "msd2", "D_r", "v2", "err"], extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
