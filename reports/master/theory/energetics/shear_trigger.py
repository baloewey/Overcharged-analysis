#!/usr/bin/env python3
"""Does the local near-wall tangential flow (shear) change before an escaper untilts? (2026-10-11). Fineflow, depleting escapes, lags -30..+2 frames.
At the escaper's wall point j (nearest wall point at the departure frame, fixed): v_t = flow along the wall tangent, depth 1-5, |along| < 4 cells, signed by the
instantaneous cavity-wide wall circulation (mean v_t over the whole wall band; avoids regression-to-mean from self-signing); also |sin tau| and cos tau of the escaper; control = same, at random wall +1/2 at the same frames."""
import os, sys, numpy as np
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
LAGS = list(range(-30, 3)); R = []; C = []
for run in sys.argv[1:]:
    J, bc, k, t, fr, ev = En.events(run, max_events=50); Q = bc["Q"]; sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Tn /= np.abs(Tn); Nn = -1j * Tn
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; D, JJ = tree.query(np.c_[X, Y])
    def vt_at(tt, j):
        ii = idx[ft[np.argmin(np.abs(ft - tt))]]; wp = W[j]; rel = (X[ii] - wp[0]) * Tn[j].real + (Y[ii] - wp[1]) * Tn[j].imag
        m = (D[ii] > 1) & (D[ii] < 5) & (np.abs(rel) < 4); loc = np.mean(rows[ii, 4][m] * Tn[j].real + rows[ii, 5][m] * Tn[j].imag) if m.any() else np.nan
        g = (D[ii] > 1) & (D[ii] < 5); Tg = Tn[JJ[ii][g]]; glob = np.mean(rows[ii, 4][g] * Tg.real + rows[ii, 5][g] * Tg.imag)
        return (loc, glob)
    def tilt(p, m):
        f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - p).T))); _, j = tree.query(p); return np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j])))))
    rng = np.random.default_rng(3)
    for e in ev:
        n = e["n"]; ch = e["chain"]; _, j = tree.query(ch[n])
        rec = {lag: (vt_at(t[n + lag], j), tilt(ch[n + lag], n + lag) if (n + lag) in ch else np.nan) for lag in LAGS if n + lag >= 0}
        R += [(l, np.sign(v[1]) * v[0], c, abs(v[1])) for l, (v, c) in rec.items()]
        f = fr[n - 20]; cand = np.flatnonzero((np.sign(f[:, 2]) == sg) & (tree.query(f[:, :2])[0] < 4))
        if len(cand):
            jc = tree.query(f[rng.choice(cand), :2])[1]; vv = {l: vt_at(t[n + l], jc) for l in LAGS if n + l >= 0}
            C += [(l, np.sign(v[1]) * v[0]) for l, v in vv.items()]
R = np.array(R); C = np.array(C); dt = 2
L = ["lag(t.u.)  v_t(escaper site, signed by global wall circulation)  cos_tau   v_t(control site)   |global wall circulation|"]
for l in LAGS:
    a = R[R[:, 0] == l]; b = C[C[:, 0] == l]
    L.append(f"{l * dt:+5d}  {np.nanmean(a[:, 1]):+.3f} ± {np.nanstd(a[:, 1]) / np.sqrt(len(a)):.3f}   {np.nanmean(a[:, 2]):+.2f}   {np.nanmean(b[:, 1]):+.3f} ± {np.nanstd(b[:, 1]) / np.sqrt(max(len(b), 1)):.3f}   {np.nanmean(a[:, 3]):.3f}")
print("\n".join(L)); open(f"shear_trigger{os.environ.get('TAG', '')}.txt", "w").write("\n".join(L) + "\n")
