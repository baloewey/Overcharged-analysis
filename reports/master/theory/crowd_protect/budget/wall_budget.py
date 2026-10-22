#!/usr/bin/env python3
"""Full charge budget of the wall band (user list, 2026-10-22): why does the steady bulk charge hold while exits slow ~x50-100 with crowding?
Both signs tracked frame to frame (2 t.u.; same-sign mutual nearest neighbours < 3.5 cells, up to 2 missed frames bridged). Band = d < 4 cells.
Per run, events per wall length per 1000 t.u.:
  P_in, P_out : +1/2 tracks crossing the band edge inward / outward;  M_in, M_out : the same for -1/2;
  C_band, A_band : pair creations / annihilations inside the band (a +1/2 and a -1/2 track born / ending in the same frame within 4 cells of each other,
                   both in the band); C_bulk, A_bulk elsewhere;
  Pb, Pd, Mb, Md : unpaired births / deaths in the band (tracker gaps, merges at the wall).
Balance of the wall +1/2 count: P_in - P_out + C_band - A_band + Pb - Pd = 0 (steady state). Second half of each run.
-> wall_budget.csv"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/crowd_protect/exchange"); import exchange_count as X  # noqa  (geometry, read_with_angle)
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__)); RL, GAP = 3.5, 2
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr); dt = (t[1] - t[0]) if nf > 1 else 2.0
        act = {+1: {}, -1: {}}; nid = 0; last = {}; first = {}; prevd = {}
        cnt = dict(P_in=0, P_out=0, M_in=0, M_out=0, C_band=0, A_band=0, C_bulk=0, A_bulk=0, Pb=0, Pd=0, Mb=0, Md=0); births = []; deaths = []
        alive_prev = set()
        for n, f in enumerate(fr):
            alive = set()
            for sg in (+1, -1):
                p = f[np.sign(f[:, 2]) == sg] if len(f) else np.zeros((0, 4)); ids = np.full(len(p), -1)
                live = [(i, xy) for i, (xy, ln) in act[sg].items() if n - ln <= 1 + GAP]
                if live and len(p):
                    LX = np.array([xy for _, xy in live]); LI = np.array([i for i, _ in live]); dd, jj = cKDTree(LX).query(p[:, :2]); _, kk = cKDTree(p[:, :2]).query(LX)
                    ok = (dd < RL) & (kk[jj] == np.arange(len(p))); ids[ok] = LI[jj[ok]]
                new = ids < 0; ids[new] = nid + np.arange(new.sum()); nid += new.sum()
                d = dist(p[:, :2]) if len(p) else np.zeros(0)
                for i, xy, dv, nw in zip(ids, p[:, :2], d, new):
                    i = int(i)
                    if nw:
                        if n > 0: births.append((n, sg, xy, dv))
                        first[i] = n
                    else:
                        was = prevd[i] < 4; now = dv < 4
                        if was and not now: cnt["P_out" if sg > 0 else "M_out"] += 1
                        elif now and not was: cnt["P_in" if sg > 0 else "M_in"] += 1
                    prevd[i] = dv; act[sg][i] = (xy, n); last[i] = (n, sg, xy, dv); alive.add(i)
                act[sg] = {i: v for i, v in act[sg].items() if n - v[1] <= 1 + GAP}
            for i in [i for i in list(last) if last[i][0] == n - 1 - GAP and i not in alive]:          # track ended GAP frames ago
                deaths.append(last.pop(i))
        for i, v in last.items():
            if v[0] < nf - 1 - GAP: deaths.append(v)
        def pairs(ev, cP, cB, sP, sM):
            byn = {}
            for e in ev: byn.setdefault(e[0], []).append(e)
            for n, es in byn.items():
                pp = [e for e in es if e[1] > 0]; mm = [e for e in es if e[1] < 0]; used = set()
                for a in pp:
                    best = None
                    for j, b in enumerate(mm):
                        if j in used: continue
                        r = np.hypot(*(np.array(a[2]) - np.array(b[2])))
                        if r < 4 and (best is None or r < best[0]): best = (r, j)
                    if best:
                        used.add(best[1]); b = mm[best[1]]
                        cnt[cP if (a[3] < 4 and b[3] < 4) else cB] += 1
                    elif a[3] < 4: cnt[sP] += 1
                for j, b in enumerate(mm):
                    if j not in used and b[3] < 4: cnt[sM] += 1
        pairs(births, "C_band", "C_bulk", "Pb", "Mb"); pairs([(e[0] + 1, e[1], e[2], e[3]) for e in deaths], "A_band", "A_bulk", "Pd", "Md")
        T = nf * dt; o = dict(run=os.path.relpath(run, ROOT), kind=kind, Lw=Lw, T=T, **info)
        for k_, v in cnt.items(): o[k_] = 1000 * v / Lw / T
        o["balance_P"] = o["P_in"] - o["P_out"] + o["C_band"] - o["A_band"] + o["Pb"] - o["Pd"]
        return o
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    runs += sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/channel_qscan/W40_*/DONE") + glob.glob(ROOT + "/channel/W40_m6_L192_s*/DONE")
                   + glob.glob(ROOT + "/channel_fix/crowd_W40_*_sym_s*/DONE") + glob.glob(ROOT + "/channel_qcross/W40_*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "kind", "Lw", "T", "k", "m", "rho", "P_in", "P_out", "M_in", "M_out", "C_band", "A_band", "C_bulk", "A_bulk", "Pb", "Pd", "Mb", "Md", "balance_P", "err"]
    with open(os.path.join(H, "wall_budget.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r), [r["err"] for r in res if "err" in r][:3])
