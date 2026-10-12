#!/usr/bin/env python3
"""Test A (2026-10-11): does the near-wall tangential flow lead the tilt of the wall defects, or the reverse?
Runs with flowfield every 20 t.u. and defects every 2 t.u. (no imposed swirl: runs/actlayout, runs/uftrack). At 24 fixed points along the wall, every flow frame (second half):
  F(t) = local flow along the wall (v_t, depth 2-6, +-4 cells along); for the nearest wall +1/2 within 4 cells along the wall: tilt tau (tail - outward normal);
  series: |v_t|, |sin tau| (tilt magnitude), alignment a = sign(tau) * sign(v_t) * |sin tau| (tilt along the flow).
Cross-correlation of fluctuations (per site, then averaged) C_XY(lag) = <dX(t) dY(t + lag)>, lag in frames (20 t.u.): peak at lag > 0 -> X leads Y.
Pairs: (|v_t|, |sin tau|), (|v_t|, a), and (signed v_t, signed tau). Usage: lead_lag.py RUNDIR [RUNDIR ...] -> prints, and lead_lag_<tag>.txt"""
import glob, gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
LAGS = range(-30, 31, 3)   # in flow frames; with 2 t.u. frames: -60 .. +60 t.u. in steps of 6 t.u.
def series(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"; op = gzip.open if fn.endswith(".gz") else open
    rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit())); ft = np.unique(rows[:, 0]); ft = ft[ft > ft.max() / 2]
    x = rows[:, 1] + 0.5 - Q[0]; y = rows[:, 2] + 0.5 - Q[1]; d, j = tree.query(np.column_stack([x, y])); nn = Nn[j]; vt_all = -rows[:, 4] * nn.imag + rows[:, 5] * nn.real
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); t = np.array(t); sites = np.linspace(0, len(W), 25)[:-1].astype(int)
    S = {q: [] for q in sites}
    for tt in ft:
        m = rows[:, 0] == tt; k = np.argmin(np.abs(t - tt))
        f = fr[k] if abs(t[k] - tt) <= 2 else np.zeros((0, 4))
        if len(f): dd, jj = tree.query(f[:, :2]); wl = np.flatnonzero((f[:, 2] > 0) & (dd < 4))
        else: wl = []
        for q in sites:
            sel = m & (d > 2) & (d < 6) & (np.abs(((s[j] - s[q] + P / 2) % P) - P / 2) < 4)
            vt = vt_all[sel].mean() if sel.any() else np.nan; tau = np.nan
            if len(wl):
                ds = np.abs(((s[jj[wl]] - s[q] + P / 2) % P) - P / 2); kk = np.argmin(ds)
                if ds[kk] < 4: i = wl[kk]; tau = np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[jj[i]]))))
            S[q].append((vt, tau))
    return {q: np.array(v) for q, v in S.items()}
def xcorr(a, b):
    out = []
    for lag in LAGS:
        if lag >= 0: x, y = a[:len(a) - lag], b[lag:]
        else: x, y = a[-lag:], b[:len(b) + lag]
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 10: out.append(np.nan); continue
        x, y = x[m] - x[m].mean(), y[m] - y[m].mean(); out.append(np.mean(x * y) / max(np.std(x) * np.std(y), 1e-12))
    return np.array(out)
if __name__ == "__main__":
    runs = sys.argv[1:]; tag = os.path.basename(os.path.dirname(os.path.abspath(runs[0]))) if runs else "x"; L = []
    for run in runs:
        try: S = series(run)
        except Exception as ex: L.append(f"{run}: {ex}"); continue
        C1, C2, C3 = [], [], []
        for q, v in S.items():
            vt, tau = v[:, 0], v[:, 1]; st = np.abs(np.sin(tau)); al = np.sign(tau) * np.sign(vt) * st
            C1.append(xcorr(np.abs(vt), st)); C2.append(xcorr(np.abs(vt), al)); C3.append(xcorr(vt, tau))
        C1, C2, C3 = (np.nanmean(c, 0) for c in (C1, C2, C3))
        L.append(f"{os.path.relpath(run, ROOT)}: lags {list(LAGS)} x (flow frame spacing) (lag > 0: flow leads tilt)")
        for nm, c in (("|v_t| vs |sin tau|", C1), ("|v_t| vs alignment", C2), ("v_t vs tau (signed)", C3)):
            L.append(f"   {nm:20s}: " + " ".join(f"{x:+.2f}" for x in c) + f"   peak at lag {list(LAGS)[int(np.nanargmax(np.abs(c)))]:+d}")
    print("\n".join(L)); open(f"lead_lag_{tag}.txt", "w").write("\n".join(L) + "\n")
