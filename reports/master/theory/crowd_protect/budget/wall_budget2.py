#!/usr/bin/env python3
"""Full charge budget of the wall band, both signs, with the validated tracker (NEW v2: crowd_protect/tracker_validation/trackers.track_new, gate 9,
2-frame gaps; user, 2026-10-22: 'recheck the wall budget with the new tracker'). Definitions as wall_budget.py: band d < 4; per wall length per 1000 t.u.:
P_in/P_out, M_in/M_out = band-edge crossings of +1/2 / -1/2 tracks (consecutive positions, gaps included); C_band/A_band = pair creations/annihilations
inside the band (a +1/2 and a -1/2 track born/ending in the same frame within 5 cells, both in the band); C_bulk/A_bulk elsewhere; Pb/Pd/Mb/Md = unpaired
track births/deaths in the band. Balance of the band's +1/2 count: P_in - P_out + C_band - A_band + Pb - Pd = 0. Exact check: net charge carried across
d = 4 by tracks + straddling pairs vs the change of the charge beyond d = 4 (rms over the run). -> wall_budget2.csv"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
ROOT = X.ROOT; H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr); dt = (t[1] - t[0]) if nf > 1 else 2.0
        L = None
        if kind == "channel": import json; L = json.load(open(run + "/input.json"))["domain"][0]
        T = track_new(fr, 9.0, L=L); c = dict(P_in=0, P_out=0, M_in=0, M_out=0, C_band=0, A_band=0, C_bulk=0, A_bulk=0, Pb=0, Pd=0, Mb=0, Md=0)
        cum = np.zeros(nf); births, deaths = [], []
        for tr in T:
            d = dist(np.array(tr["xy"])); inb = d < 4; s = "P" if tr["sg"] > 0 else "M"
            for i in np.flatnonzero(inb[1:] != inb[:-1]):
                c[s + ("_in" if inb[i + 1] else "_out")] += 1; cum[tr["n"][i + 1]] += 0.5 * tr["sg"] * (-1 if inb[i + 1] else 1)   # charge beyond d = 4
            if tr["n"][0] > 0: births.append((tr["n"][0], tr["sg"], np.array(tr["xy"][0]), d[0]))
            if tr["n"][-1] < nf - 3: deaths.append((tr["n"][-1] + 1, tr["sg"], np.array(tr["xy"][-1]), d[-1]))
        for ev, cB, cO, sP, sM, sgn in ((births, "C_band", "C_bulk", "Pb", "Mb", 1), (deaths, "A_band", "A_bulk", "Pd", "Md", -1)):
            byn = {}
            for e in ev: byn.setdefault(e[0], []).append(e)
            for n, es in byn.items():
                pp = [e for e in es if e[1] > 0]; mm = [e for e in es if e[1] < 0]; used = set()
                for a in pp:
                    best = None
                    for j, b in enumerate(mm):
                        if j in used: continue
                        dx = a[2][0] - b[2][0]
                        if L: dx = (dx + L / 2) % L - L / 2
                        r = np.hypot(dx, a[2][1] - b[2][1])
                        if r < 5 and (best is None or r < best[0]): best = (r, j)
                    if best:
                        used.add(best[1]); b = mm[best[1]]; inA, inB = a[3] < 4, b[3] < 4
                        c[cB if (inA and inB) else cO] += 1
                        if inA != inB: cum[min(n, nf - 1)] += sgn * (0.5 if not inA else -0.5)          # straddling pair moves charge across d = 4
                    elif a[3] < 4: c[sP] += 1
                for j, b in enumerate(mm):
                    if j not in used and b[3] < 4: c[sM] += 1
        ex = np.array([0.5 * np.sign(f[dist(f[:, :2]) > 4, 2]).sum() if len(f) else np.nan for f in fr]); ex -= np.nanmean(ex[:10])
        Tt = nf * dt; o = dict(run=os.path.relpath(run, ROOT), kind=kind, Lw=Lw, T=Tt, **info, closure_rms=float(np.sqrt(np.nanmean((np.cumsum(cum) - ex) ** 2))))
        for k_, v in c.items(): o[k_] = 1000 * v / Lw / Tt
        o["balance_P"] = o["P_in"] - o["P_out"] + o["C_band"] - o["A_band"] + o["Pb"] - o["Pd"]
        return o
    except Exception as ex_:
        return dict(run=os.path.relpath(run, ROOT), err=str(ex_))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    runs += sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/channel_qscan/W40_*/DONE") + glob.glob(ROOT + "/channel/W40_m6_L192_s*/DONE") + glob.glob(ROOT + "/channel_fix/crowd_W40_*_sym_s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "kind", "Lw", "T", "k", "m", "rho", "P_in", "P_out", "M_in", "M_out", "C_band", "A_band", "C_bulk", "A_bulk", "Pb", "Pd", "Mb", "Md", "balance_P", "closure_rms", "err"]
    with open(os.path.join(H, "wall_budget2.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r), [r["err"] for r in res if "err" in r][:2])
