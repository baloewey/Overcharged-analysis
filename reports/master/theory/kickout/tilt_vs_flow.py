#!/usr/bin/env python3
"""Empirical tilt-flow relation for wall +1/2: cos(tau) (and the signed tilt) vs the local flow along the wall at the site (|v_t| and signed v_t, depth 2-6, +-4 cells
along the wall, flowfield every 20 t.u.), in the imposed-swirl runs (passive and active). Also the autocorrelation time of the local v_t at fixed wall points.
Input to the combined escape model. -> tilt_vs_flow.txt"""
import glob, gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); VB = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 1.0])
def analyse(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"; op = gzip.open if fn.endswith(".gz") else open
    rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit())); ft = np.unique(rows[:, 0])
    x = rows[:, 1] + 0.5 - Q[0]; y = rows[:, 2] + 0.5 - Q[1]; d, j = tree.query(np.column_stack([x, y])); nn = Nn[j]
    vt_all = -rows[:, 4] * nn.imag + rows[:, 5] * nn.real
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); t = np.array(t); out = []; tseries = {}
    Tend = J["simSteps"] * 0.1; tf = J.get("swirlFlipTime", 0) or 1e9
    for tt in ft:
        if tt < Tend / 2 or abs(tt - tf) < 300: continue
        k = np.argmin(np.abs(t - tt))
        if abs(t[k] - tt) > 2 or not len(fr[k]): continue
        f = fr[k]; dd, jj = tree.query(f[:, :2]); w = (f[:, 2] > 0) & (dd < 4); m = rows[:, 0] == tt
        for i in np.flatnonzero(w):
            sel = m & (d > 2) & (d < 6) & (np.abs(((s[j] - s[jj[i]] + P / 2) % P) - P / 2) < 4)
            if not sel.any(): continue
            vt = vt_all[sel].mean(); tau = np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[jj[i]]))))
            out.append((abs(vt), np.cos(tau), np.sign(vt) * np.sign(tau)))
        for q in range(0, len(W), max(1, len(W) // 24)):                      # fixed wall points for the autocorrelation
            sel = m & (d > 2) & (d < 6) & (np.abs(((s[j] - s[q] + P / 2) % P) - P / 2) < 4)
            if sel.any(): tseries.setdefault(q, []).append(vt_all[sel].mean())
    return np.array(out), tseries, ft[1] - ft[0]
L = []
for lab, pat in (("passive swirl", "swirl_ellipse/passive_*"), ("active swirl", "swirl_ellipse*/active_*")):
    O = []; acs = []; dtf = None
    for run in sorted(glob.glob(os.path.join(ROOT, pat))):
        try: o, ts, dtf = analyse(run); O.append(o)
        except Exception as ex: print(run, ex); continue
        for v in ts.values():
            v = np.array(v) - np.mean(v)
            if len(v) > 10 and v.var() > 0: acs.append([np.mean(v[l:] * v[:len(v) - l]) / v.var() for l in range(6)])
    O = np.concatenate(O); L.append(f"{lab}: {len(O)} wall-defect samples; cos(tau) by local |v_t| bin:")
    for a, b in zip(VB[:-1], VB[1:]):
        m = (O[:, 0] >= a) & (O[:, 0] < b)
        if m.sum() > 20: L.append(f"   |v_t| {a:.2f}-{b:.2f}: cos tau {O[m, 1].mean():+.3f}  tilt along the flow {O[m, 2].mean():+.2f}  (n {m.sum()})")
    if acs: L.append(f"   local v_t autocorrelation at lags 0..5 x {dtf:.0f} t.u.: " + " ".join(f"{x:.2f}" for x in np.mean(acs, 0)))
print("\n".join(L)); open(os.path.join(H, "tilt_vs_flow.txt"), "w").write("\n".join(L) + "\n")
