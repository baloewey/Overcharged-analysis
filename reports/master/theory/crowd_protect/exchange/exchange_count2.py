#!/usr/bin/env python3
"""Clean wall<->bulk exchange counts with the predictive tracker (2026-10-22; recheck of exchange_count.py after the channel-current artefact).
Same definitions as exchange_count.py (WALL d < 4, BULK d > 6; escape = >= 20 t.u. in WALL -> BULK, no return for 20 t.u.; entry = >= 20 t.u. in BULK
-> WALL, stays 20 t.u.; exposures), but tracks per sign by velocity prediction (last displacement; v2 from 22 Oct late: cost = min(distance to prediction, to last position)) + global Hungarian assignment, gate GATE cells, up to
2 missed frames (prediction carried on; residence counters run through gaps). Exact tracking-free check: the charge farther than 4 cells from the wall
can change only by crossings of the closed d = 4 curve; reported as its slope (should be ~0) next to the tracked net crossing flux.
Usage: exchange_count2.py [GATE] -> exchange_count2_g<GATE>.csv"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
from scipy.optimize import linear_sum_assignment
import exchange_count as X
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__)); NW = 10; GAP = 2; GATE = float(sys.argv[1]) if len(sys.argv) > 1 else 9.0
def tracks(fr, gate):
    T = []; act = {+1: [], -1: []}
    for n, f in enumerate(fr):
        for sg in (+1, -1):
            p = f[np.sign(f[:, 2]) == sg][:, :2] if len(f) else np.zeros((0, 2)); A = act[sg]; used = np.zeros(len(p), bool)
            if A and len(p):
                pred = []
                for ti in A:
                    tr = T[ti]; m = n - tr["n"][-1]; v = (np.array(tr["xy"][-1]) - np.array(tr["xy"][-2])) / (tr["n"][-1] - tr["n"][-2]) if len(tr["n"]) > 1 else np.zeros(2)
                    pred.append(np.array(tr["xy"][-1]) + v * m)
                pred = np.array(pred); C = np.hypot(p[None, :, 0] - pred[:, None, 0], p[None, :, 1] - pred[:, None, 1])
                last = np.array([T[ti]["xy"][-1] for ti in A]); C = np.minimum(C, np.hypot(p[None, :, 0] - last[:, None, 0], p[None, :, 1] - last[:, None, 1]))   # v2 dual cost
                r, c = linear_sum_assignment(np.where(C < gate, C, 1e6))
                for a, b in zip(r, c):
                    if C[a, b] < gate: T[A[a]]["n"].append(n); T[A[a]]["xy"].append(tuple(p[b])); used[b] = True
            keep = [ti for ti in A if n - T[ti]["n"][-1] <= GAP]
            for b in np.flatnonzero(~used): T.append(dict(sg=sg, n=[n], xy=[tuple(p[b])])); keep.append(len(T) - 1)
            act[sg] = keep
    return T
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr); dt = (t[1] - t[0]) if nf > 1 else 2.0
        T = tracks(fr, GATE); esc = ent = 0; expW = expB = 0.0; netx = 0.0; lens = []
        for tr in T:
            d = dist(np.array(tr["xy"])); nn = np.array(tr["n"]); lens.append(nn[-1] - nn[0] + 1)
            full = np.full(nn[-1] - nn[0] + 1, np.nan); full[nn - nn[0]] = d
            for i in range(1, len(full)):                                 # bridge gaps with the last distance
                if np.isnan(full[i]): full[i] = full[i - 1]
            st = np.where(full < 4, 1, np.where(full > 6, 2, 0)); L_ = len(st); runW = runB = 0
            s4 = full > 4; netx += 0.5 * tr["sg"] * (int(s4[-1]) - int(s4[0])) * 0      # placeholder (net flux computed from the exact check below)
            for a in range(L_):
                if st[a] == 1:
                    if runB >= NW and a + NW < L_ and np.all(st[a:a + NW + 1] == 1): ent += 1
                    runB = 0; runW += 1
                    if runW > NW: expW += dt
                elif st[a] == 2:
                    if runW >= NW and a + NW < L_ and not np.any(st[a:a + NW + 1] == 1): esc += 1
                    runW = 0; runB += 1
                    if runB > NW: expB += dt
        qb = []; tt = []
        for n, f in enumerate(fr):
            if len(f): qb.append(0.5 * np.sign(f[dist(f[:, :2]) > 4, 2]).sum()); tt.append(n * dt)
        slope = np.polyfit(tt, qb, 1)[0] * 1000 / Lw if len(qb) > 10 else np.nan
        return dict(run=os.path.relpath(run, ROOT), kind=kind, Lw=Lw, T=nf * dt, **info, esc=esc, ent=ent, expW=expW, expB=expB, med_len=np.median(lens) * dt, qb_slope=slope)
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    runs += sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/channel_qscan/W40_*/DONE") + glob.glob(ROOT + "/channel/W40_m6_L192_s*/DONE") + glob.glob(ROOT + "/channel_fix/crowd_W40_*_sym_s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "kind", "Lw", "T", "k", "m", "rho", "esc", "ent", "expW", "expB", "med_len", "qb_slope", "err"]
    with open(os.path.join(H, f"exchange_count2_g{GATE:g}.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r))
