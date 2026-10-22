#!/usr/bin/env python3
"""Clean wall<->bulk exchange counts (user, 2026-10-22: 'do the defect counting'). Every +1/2 tracked frame to frame (2 t.u.; mutual nearest neighbours
< 3.5 cells, up to 2 missed frames bridged). Distance d to the nearest wall (cavity: wall polyline; channel: min(y - y_bottom, y_top - y)). States: WALL d < 4, BULK d > 6.
  ESCAPE: a track continuously in WALL for >= 20 t.u. reaches BULK and does not return to WALL for the next 20 t.u. (track must persist).
  ENTRY : a track continuously in BULK for >= 20 t.u. reaches WALL and stays in WALL for the next 20 t.u.
  Exposures: wall-resident time (track in WALL with >= 20 t.u. history) and bulk-resident time (in BULK with >= 20 t.u. history), defect x t.u.
  Also: tracks born / ending inside WALL (pair creation, annihilation, tracker gaps), mean wall and bulk +1/2 counts, bulk charge (+-1/2 with d > 6).
Second half of each run. Shape sweep (675) + channels (channel_qscan, channel W40_m6, channel_fix crowd sym, channel_qcross when done). -> exchange_count.csv"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); NW = 10; RL = 3.5; GAP = 2    # 10 frames = 20 t.u.; link radius (cells); missed frames bridged
def n0_gap(e): return 0
def geometry(run):
    J = json.load(open(run + "/input.json")); BC = J["BC"]
    if J["domain"][1] < 80 and len(BC) >= 4 and BC[2].get("customAnch") is not None and BC[2].get("aInv", [0, 0, 0])[:2] == [0, 1]:
        y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; L = J["domain"][0]; m = abs(BC[2]["customAnch"])
        return "channel", np.array([0.0, 0.0]), (lambda xy: np.minimum(xy[:, 1] - y0, y1 - xy[:, 1])), 2.0 * L, dict(m=m, L=L, W=y1 - y0, rho=2 * m / L)
    bc = BC[0]; W, s, P = wall_polyline(bc); tree = cKDTree(W)
    return "cavity", np.array(bc["Q"][:2], float), (lambda xy: tree.query(xy)[0]), P, dict(k=bc["customAnch"], rho=2 * bc["customAnch"] / P)
def one(run):
    try:
        kind, Q, dist, Lw, info = geometry(run); t, fr = read_with_angle(run, Q, 0.5); nf = len(fr); dt = (t[1] - t[0]) if nf > 1 else 2.0
        nid = 0; D = {}; act = {}                                       # act: id -> (last xy, last frame)
        nW, nB, qb = [], [], []
        for n, f in enumerate(fr):
            if len(f):
                dall = dist(f[:, :2]); qb.append(0.5 * np.sign(f[dall > 6, 2]).sum())
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)); ids = np.full(len(p), -1)
            live = [(i_, xy) for i_, (xy, ln) in act.items() if n - ln <= 1 + GAP]
            if live and len(p):
                LX = np.array([xy for _, xy in live]); LI = np.array([i_ for i_, _ in live])
                dd, jj = cKDTree(LX).query(p[:, :2]); _, kk = cKDTree(p[:, :2]).query(LX); ok = (dd < RL) & (kk[jj] == np.arange(len(p))); ids[ok] = LI[jj[ok]]
            new = ids < 0; ids[new] = nid + np.arange(new.sum()); nid += new.sum()
            for i_, xy in zip(ids, p[:, :2]): act[int(i_)] = (xy, n)
            act = {i_: v for i_, v in act.items() if n - v[1] <= 1 + GAP}
            d = dist(p[:, :2]) if len(p) else np.zeros(0); nW.append((d < 4).sum()); nB.append((d > 6).sum())
            for i_, dv in zip(ids, d):
                e = D.setdefault(int(i_), [n, []])
                while n0_gap(e) < n - e[0] - len(e[1]): e[1].append(e[1][-1])          # bridge missed frames with the last distance
                e[1].append(dv)
        esc = ent = born_w = end_w = 0; expW = expB = 0.0
        for i_, (n0, ds) in D.items():
            ds = np.array(ds); st = np.where(ds < 4, 1, np.where(ds > 6, 2, 0)); nn = len(st)
            if st[0] == 1 and n0 > 0: born_w += 1
            if st[-1] == 1 and n0 + nn < nf: end_w += 1
            runW = runB = 0                                               # residence counters survive the 4-6 cell band
            for a in range(nn):
                if st[a] == 1:
                    if runB >= NW and a + NW < nn and np.all(st[a:a + NW + 1] == 1): ent += 1
                    runB = 0; runW += 1
                    if runW > NW: expW += dt
                elif st[a] == 2:
                    if runW >= NW and a + NW < nn and not np.any(st[a:a + NW + 1] == 1): esc += 1
                    runW = 0; runB += 1
                    if runB > NW: expB += dt
        T = nf * dt
        return dict(run=os.path.relpath(run, ROOT), kind=kind, Lw=Lw, T=T, **info, esc=esc, ent=ent, expW=expW, expB=expB, born_w=born_w, end_w=end_w,
                    nW=np.mean(nW), nB=np.mean(nB), qb=np.mean(qb) if qb else np.nan)
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    runs += sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/channel_qscan/W40_*/DONE") + glob.glob(ROOT + "/channel/W40_m6_L192_s*/DONE")
                   + glob.glob(ROOT + "/channel_fix/crowd_W40_*_sym_s*/DONE") + glob.glob(ROOT + "/channel_qcross/W40_*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "kind", "Lw", "T", "k", "m", "L", "W", "rho", "esc", "ent", "expW", "expB", "born_w", "end_w", "nW", "nB", "qb", "err"]
    with open(os.path.join(H, "exchange_count.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r), [r["err"] for r in res if "err" in r][:3])
