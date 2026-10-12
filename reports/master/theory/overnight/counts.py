#!/usr/bin/env python3
"""Overnight analysis of the imposed-swirl batches (runs/swirl_neg, swirl_ellipse, swirl_threshold): per run and time window, wall +1/2 and -1/2 counts,
bulk charge, circulation (mean v_theta, normalised by v_rms), two-vortex fraction (halves along x), tilt coherence of the wall +1/2, escapes (wall +1/2
followed 20 t.u. ending deeper than 6 cells). Windows: start (0-200), before the flip, after the flip (or the second half without a flip).
Usage: counts.py BATCHDIR > BATCHDIR_summary.txt"""
import glob, gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.path as mpath
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
batch = sys.argv[1]
for run in sorted(glob.glob(os.path.join(ROOT, batch, "*"))):
    if not os.path.isdir(run) or not os.path.exists(run + "/DONE"): continue
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn); path = mpath.Path(W)
    tf = J.get("swirlFlipTime", 0.0); Tend = J["simSteps"] * 0.1
    wins = [("start", 0, 200)] + ([("before flip", tf - 500, tf), ("after flip", Tend - 500, Tend)] if tf > 0 else [("second half", Tend / 2, Tend)])
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); t = np.array(t)
    fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"
    op = gzip.open if fn.endswith(".gz") else open
    rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit())) if os.path.exists(fn) else None
    line = f"{os.path.basename(run):28s} (k {k}, swirl {J.get('swirlAcc', 0)}, tumble {J['species'][0].get('tumble', 2)}):"
    for lab, a, b in wins:
        sel = [i for i in range(len(t)) if a <= t[i] < b]
        wp, wm, qb, coh = [], [], [], []
        for i in sel[::5]:
            f = fr[i]; d, j = tree.query(f[:, :2]) if len(f) else (np.array([]), np.array([], int))
            pl = f[:, 2] > 0 if len(f) else np.array([], bool)
            wp.append(((d < 4) & pl).sum()); wm.append(((d < 4) & ~pl).sum()); qb.append(0.5 * np.sign(f[d > 6, 2]).sum() if len(f) else 0)
            w = (d < 4) & pl
            if w.any(): coh.append(np.mean(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))))
        circ = two = np.nan
        if rows is not None:
            ts = np.unique(rows[:, 0]); ts = ts[(ts > a) & (ts <= b)]; cs, tw = [], []
            for tt in ts:
                m = rows[:, 0] == tt; x = rows[m, 1] + 0.5 - Q[0]; y = rows[m, 2] + 0.5 - Q[1]; vx, vy = rows[m, 4], rows[m, 5]
                ins = path.contains_points(np.column_stack([x, y])); r = np.hypot(x, y); vrms = np.sqrt(np.mean((vx ** 2 + vy ** 2)[ins]))
                cs.append(np.mean(((x * vy - y * vx) / np.maximum(r, 0.5))[ins]) / max(vrms, 1e-9))
                def hr(sl):
                    cx, cy = x[sl].mean(), y[sl].mean(); dx, dy = x[sl] - cx, y[sl] - cy
                    return np.mean((dx * vy[sl] - dy * vx[sl]) / np.maximum(np.hypot(dx, dy), 0.5))
                tw.append(hr(ins & (x < 0)) * hr(ins & (x > 0)) < 0)
            circ = np.mean(cs) if cs else np.nan; two = np.mean(tw) if tw else np.nan
        line += (f" | {lab}: wall +{np.mean(wp):.1f} (mode {np.bincount(wp).argmax() if wp else -1}) -{np.mean(wm):.1f}, Q_bulk {np.mean(qb):+.2f}, "
                 f"circ/vrms {circ:+.2f}, P(two) {two:.2f}, tilt coh {np.mean(np.abs(coh)) if coh else np.nan:.2f}")
    print(line, flush=True)
