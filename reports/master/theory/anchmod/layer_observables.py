#!/usr/bin/env python3
"""Geometry vs topology inside the circle (user's question, 2026-10-07): does moving the pinned sites (modulated winding k phi + A sin(m phi))
change the layer? Per run (k = 6, R = 31.83, active), second half of the run:
  topology   bulk charge Q_bulk = sum of charges deeper than 6 cells (mean, and fraction of frames with Q_bulk = 1); wall +1/2 count (mode,
             fraction of frames with 2k - 2 = 10); vacancy fraction
  layer      |tilt| (circular mean of the wall tails vs normal), tilt coherence, reversal rate (sign flips of the mean tilt per 1000 t.u.,
             the tilt sign follows the circulation)
  survival   escape rate of wall +1/2: transitions wall (depth < 4) -> bulk (depth > 6) per wall defect per 1000 t.u. (nearest-neighbour
             tracking between consecutive frames, dt = 2)
  flow       (if flowfield present) |circulation|, fraction of intervals with two vortices (halves about the x and the y axis counter-rotating)
Groups: baseline A = 0 (Rscan_R31.83_seeds s01-s03, pinning_ic_active, anchmod_sweep A0), anchmod_test, anchmod_sweep.
Writes layer_observables.txt."""
import glob, gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
def analyse(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = int(bc["customAnch"])
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); t = np.array(t)
    nw, qb, tl, esc, wall_time = [], [], [], 0, 0.0
    prev = None
    for i, f in enumerate(fr):
        dep = R - np.hypot(f[:, 0], f[:, 1]); pl = f[:, 2] > 0
        w = pl & (dep < 4); nw.append(w.sum()); qb.append(0.5 * np.sum(np.sign(f[dep > 6, 2])))
        if w.any(): tl.append(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))))
        if prev is not None:
            pw = prev[(prev[:, 2] > 0) & (R - np.hypot(prev[:, 0], prev[:, 1]) < 4)]
            cur = f[pl]
            if len(pw) and len(cur):
                d, j = cKDTree(cur[:, :2]).query(pw[:, :2])
                ok = d < 3; esc += np.sum(dep[pl][j[ok]] > 6); wall_time += len(pw) * (t[i] - t[i - 1])
        prev = f
    nw = np.array(nw); qb = np.array(qb); tl = np.array(tl); ta = np.angle(tl)
    sgn = np.sign(ta); flips = np.sum(np.diff(sgn[sgn != 0]) != 0); T = t[-1] - t[0]
    out = dict(k=k, A=bc.get("anchAmp", 0.0), m=bc.get("anchMode", 0), Qb=qb.mean(), Qb1=np.mean(np.abs(qb - 1) < 0.01), nw_mode=np.bincount(nw).argmax(),
               f10=np.mean(nw == 2 * k - 2), tilt=np.degrees(np.mean(np.abs(ta))), coh=np.mean(np.abs(tl)), flips=flips / T * 1000,
               esc=esc / max(wall_time, 1) * 1000)
    ff = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else (run + "/flowfield.dat" if os.path.exists(run + "/flowfield.dat") else None)
    if ff:
        op = gzip.open if ff.endswith(".gz") else open; Q = bc["Q"]
        rows = np.loadtxt((l for l in op(ff, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]; circ, two = [], []
        for tt in ts:
            r_ = rows[rows[:, 0] == tt]; x = r_[:, 1] + 0.5 - Q[0]; y = r_[:, 2] + 0.5 - Q[1]; rr = np.hypot(x, y); m = rr < R - 1
            vt = (x * r_[:, 5] - y * r_[:, 4]) / np.maximum(rr, 0.5); circ.append(np.mean(vt[m]))
            def halfrot(sel):
                cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
                return np.mean((dx * r_[sel, 5] - dy * r_[sel, 4]) / np.maximum(np.hypot(dx, dy), 0.5))
            tw = [halfrot(m & (x < 0)) * halfrot(m & (x > 0)) < 0, halfrot(m & (y < 0)) * halfrot(m & (y > 0)) < 0]
            two.append(any(tw))
        out.update(circ=np.mean(np.abs(circ)), two=np.mean(two))
    return out
groups = [("baseline A=0", [ROOT + f"/cavity/Rscan_R31.83_seeds/{s}/k06" for s in ("s01", "s02", "s03")] + sorted(glob.glob(ROOT + "/pinning_ic_active/*")) + sorted(glob.glob(ROOT + "/anchmod_sweep/A0.0_*"))),
          ("anchmod_test (no flow output)", sorted(glob.glob(ROOT + "/anchmod_test/active_*"))),
          ("anchmod_sweep", sorted(d for d in glob.glob(ROOT + "/anchmod_sweep/A*") if not os.path.basename(d).startswith("A0.0")))]
hdr = f"{'run':34s} {'A':>4s} {'m':>2s} | {'Q_bulk':>6s} {'P(Q=1)':>6s} {'n_wall':>6s} {'P(10)':>5s} | {'|tilt|':>6s} {'coh':>5s} {'rev/1000':>8s} | {'escape/1000':>11s} | {'|circ|':>6s} {'P(two)':>6s}"
for g, runs in groups:
    say(f"\n== {g}"); say(hdr)
    for run in runs:
        if not os.path.exists(run + "/DONE"): say(f"{os.path.basename(run):34s} not finished"); continue
        o = analyse(run)
        say(f"{run.split('runs/')[1][:34]:34s} {o['A']:4.1f} {o['m']:2.0f} | {o['Qb']:6.2f} {o['Qb1']:6.2f} {o['nw_mode']:6d} {o['f10']:5.2f} | {o['tilt']:6.1f} {o['coh']:5.2f} {o['flips']:8.2f} | "
            f"{o['esc']:11.3f} | {o.get('circ', np.nan):6.3f} {o.get('two', np.nan):6.2f}")
open(os.path.join(H, "layer_observables.txt"), "w").write("\n".join(L) + "\n")
