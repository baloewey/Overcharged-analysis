#!/usr/bin/env python3
"""Do wall escapes cascade? (user's question 2026-10-09: does an escape at the tips destabilise the layer and help the next escape?)
Shape sweep, 675 runs, second half, frames every 10 t.u. Escape event = wall +1/2 (depth < 4) whose nearest-neighbour chain (< 6 cells per 10 t.u.) is deeper
than 6 cells 20 t.u. later; repeated detections of one departure (same site within 4 cells of arclength, within 20 t.u.) are merged (latest start kept).
(1) temporal clustering: pair correlation g(tau) = observed later-event pairs at lag tau / expected for a Poisson process at the run's rate (edge-corrected),
    pooled per class; avalanche -> excess at short lags above the slow-modulation level at 500-2000 t.u.
(2) spatial: arclength separation (in site spacings P/(2k-2)) of event pairs with lag < 100 t.u., relative to pairs at lag 500-2000 t.u. in the same run
    (null with the same site preferences); elastic cascade -> excess at 1 spacing.
(3) hazard: probability that a wall defect escapes within 20 t.u. vs the instantaneous vacancy count n_v = (2k-2) - n_wall; pooled and within-run
    (n_v - run mean, p / run mean p). Runaway -> rising; restoring -> falling.
Writes avalanche.txt (+ avalanche_events.pkl)."""
import csv, glob, json, os, pickle, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max()
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t)[::5]; fr = fr[::5]
        PL, DEP, SS = [], [], []
        for f in fr:
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)); d, i = tree.query(p[:, :2]) if len(p) else (np.zeros(0), np.zeros(0, int))
            PL.append(p); DEP.append(d); SS.append(s[i] if len(p) else np.zeros(0))
        raw = []; hz = []                                   # raw events (t, s, |x|/A); hazard samples (n_v, n_wall, n_escaping)
        for a in range(len(PL) - 2):
            w = DEP[a] < 4; nwall = int(w.sum()); nesc = 0
            if nwall:
                pos = PL[a][w, :2]; s0 = SS[a][w]; x0 = np.abs(pos[:, 0]) / A; alive = np.ones(nwall, bool)
                for b in (1, 2):
                    if not len(PL[a + b]): alive[:] = False; break
                    dd, jj = cKDTree(PL[a + b][:, :2]).query(pos); alive &= dd < 6; pos = PL[a + b][jj, :2]; dfin = DEP[a + b][jj]
                if len(PL[a + 2]):
                    e = alive & (dfin > 6); nesc = int(e.sum())
                    raw += [(t[a], ss, xx) for ss, xx in zip(s0[e], x0[e])]
            hz.append((2 * k - 2 - nwall, nwall, nesc))
        ev = []                                              # merge repeated detections
        for te, se, xe in sorted(raw):
            if ev and any(te - u[0] <= 20 and min(abs(se - u[1]), P - abs(se - u[1])) < 4 for u in ev[-6:]):
                for m in range(len(ev) - 1, max(len(ev) - 7, -1), -1):
                    if te - ev[m][0] <= 20 and min(abs(se - ev[m][1]), P - abs(se - ev[m][1])) < 4: ev[m] = (te, se, xe); break
            else: ev.append((te, se, xe))
        return dict(rel=rel, k=k, P=P, T0=t[0], T1=t[-3], ev=np.array(ev).reshape(-1, 3), hz=np.array(hz))
    except Exception as ex:
        return dict(rel=rel, error=str(ex))
LAGS = [(10, 30), (30, 60), (60, 100), (100, 200), (200, 500), (500, 1000), (1000, 2000)]
def clustering(RS):
    obs = np.zeros(len(LAGS)); exp = np.zeros(len(LAGS))
    for r in RS:
        tt = r["ev"][:, 0]; n = len(tt); T = r["T1"] - r["T0"]
        if n < 2: continue
        rate = n / T
        for b, (l1, l2) in enumerate(LAGS):
            dt = tt[None, :] - tt[:, None]; obs[b] += ((dt >= l1) & (dt < l2)).sum()
            exp[b] += rate * np.sum(np.clip(np.minimum(tt + l2, r["T1"]) - np.minimum(tt + l1, r["T1"]), 0, None))
    return obs, exp
SB = [(0, 0.5), (0.5, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 99)]
def spatial(RS):
    near = np.zeros(len(SB)); far = np.zeros(len(SB))
    for r in RS:
        e = r["ev"]
        if len(e) < 2: continue
        sp = r["P"] / (2 * r["k"] - 2); dt = e[None, :, 0] - e[:, None, 0]; ds = np.abs(e[None, :, 1] - e[:, None, 1]); ds = np.minimum(ds, r["P"] - ds) / sp
        for b, (a1, a2) in enumerate(SB):
            m = (ds >= a1) & (ds < a2); near[b] += (m & (dt > 0) & (dt < 100)).sum(); far[b] += (m & (dt >= 500) & (dt < 2000)).sum()
    return near / max(near.sum(), 1), far / max(far.sum(), 1), near.sum()
def hazard(RS):
    pooled = {}; within = {}
    for r in RS:
        h = r["hz"]; h = h[h[:, 1] > 0]
        if len(h) < 50 or h[:, 2].sum() == 0: continue
        pm = h[:, 2].sum() / h[:, 1].sum(); nvm = np.mean(h[:, 0])
        for nv, nw, ne in h:
            pooled.setdefault(int(nv), [0, 0]); pooled[int(nv)][0] += ne; pooled[int(nv)][1] += nw
            dv = int(np.clip(np.round(nv - nvm), -3, 3)); within.setdefault(dv, [0, 0]); within[dv][0] += ne; within[dv][1] += nw * pm
    return pooled, within
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    pickle.dump(res, open(os.path.join(H, "avalanche_events.pkl"), "wb"))
    res = [r for r in res if "error" not in r]; L = []
    def say(x): print(x, flush=True); L.append(x)
    for r in res: m = man[r["rel"]]; r["el"] = float(m["aspect"]) >= 1.25; r["R"] = float(m["Req_la"])
    groups = [("circle-like (aspect <= 1.22)", lambda r: not r["el"]), ("elongated (aspect >= 1.25)", lambda r: r["el"])] + \
             [(f"{'elongated' if e else 'circle-like'} R_eq {R}", (lambda e, R: lambda r: r["el"] == e and r["R"] == R)(e, R)) for e in (False, True) for R in (1.0, 1.5, 2.0)]
    say("(1) TEMPORAL CLUSTERING: g(lag) = observed / Poisson-expected later escapes (1 = independent)")
    for name, sel in groups:
        RS = [r for r in res if sel(r)]; obs, exp = clustering(RS); nev = sum(len(r["ev"]) for r in RS)
        say(f"  {name:32s} ({len(RS)} runs, {nev} escapes): " + "  ".join(f"{l1}-{l2}: {o / e:.2f}" for (l1, l2), o, e in zip(LAGS, obs, np.maximum(exp, 1e-9))))
    say("\n(2) SPATIAL: arclength separation (site spacings) of escape pairs with lag < 100 t.u. vs pairs at lag 500-2000 (same runs); ratio > 1 = local cascade")
    for name, sel in groups[:2] + groups[2:]:
        RS = [r for r in res if sel(r)]; near, far, n = spatial(RS)
        say(f"  {name:32s} (n {int(n)} close pairs): " + "  ".join(f"{a1}-{a2 if a2 < 99 else 'inf'}: {x:.2f}/{y:.2f}={x / max(y, 1e-9):.2f}" for (a1, a2), x, y in zip(SB, near, far)))
    say("\n(3) HAZARD: P(a wall +1/2 escapes within 20 t.u.) vs instantaneous vacancies n_v; within-run: n_v - run mean -> p / run-mean p")
    for name, sel in groups:
        RS = [r for r in res if sel(r)]; pooled, within = hazard(RS)
        say(f"  {name:32s} pooled: " + "  ".join(f"nv{nv:+d}: {a / b:.4f}" for nv, (a, b) in sorted(pooled.items()) if b > 2000))
        say(f"  {'':32s} within: " + "  ".join(f"{dv:+d}: {a / b:.2f}" for dv, (a, b) in sorted(within.items()) if b > 20))
    # tips: does an escape at a tip raise the escape probability of the rest of the layer (non-tip sites) in the next 100 t.u.?
    say("\n(3b) after a TIP escape (|x|/A > 0.5): non-tip escapes in the next 100 t.u. / Poisson expectation from the run's non-tip rate (elongated)")
    for R in (1.0, 1.5, 2.0, None):
        RS = [r for r in res if r["el"] and (R is None or r["R"] == R)]; o = e_ = 0.0
        for r in RS:
            ev = r["ev"]; tip = ev[:, 2] > 0.5; nt = ev[~tip, 0]; rate = len(nt) / (r["T1"] - r["T0"])
            for te in ev[tip, 0]: o += ((nt > te) & (nt <= te + 100)).sum(); e_ += rate * max(min(te + 100, r["T1"]) - te, 0)
        say(f"  R_eq {R if R else 'all'}: {o / max(e_, 1e-9):.2f} (obs {int(o)}, exp {e_:.0f})")
    open(os.path.join(H, "avalanche.txt"), "w").write("\n".join(L) + "\n")
