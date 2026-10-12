#!/usr/bin/env python3
"""Wall-layer response to an imposed swirl and its reversal (user's suggestion, 2026-10-08; runs/swirl_test, binary bin/mpcd_swirl.out).
Azimuthal body acceleration A within 0.6R, reversed at production time t_f = 1500. k = 6, R = 31.83; passive (2k = 12 wall +1/2) and active (10).
Per run, time series: circulation C (mean v_theta, r < R - 1) and fluid at the layer depth V_L (1-2.5 cells from the wall), every 10 t.u.;
wall tilt tau (tail convention) and registry shift X = unwrap(arg <exp(i N phi)>) R_l/N (N = number of sites: 10 active, 12 passive) every 2 t.u.
Reported: plateau values before (1000-1500) and after (2500-3000) the flip for C, V_L, tau, X; time for C and tau to cross zero after the flip;
registry shift dX and its duration (10-90%); layer speed / fluid speed at the layer depth during the transit (effective advection); for the
passive runs the static tilt induced by the imposed flow. Writes swirl_response.txt and swirl_response.png."""
import glob, gzip, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
fig, ax = plt.subplots(2, 3, figsize=(14, 6.5)); TF = 1500.0
for d in sorted(glob.glob(ROOT + "/swirl_test/*")):
    name = os.path.basename(d)
    if not os.path.exists(d + "/DONE"): say(f"{name}: not finished"); continue
    J = json.load(open(d + "/input.json")); bc = J["BC"][0]; R = bc["R"]; Q = bc["Q"]; k = int(bc["customAnch"]); active = name.startswith("active")
    N = 2 * k - 2 if active else 2 * k; Rl = R - 1.6; row = 0 if active else 1
    fn = d + "/flowfield.dat.gz"; A = np.array([l.split() for l in gzip.open(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()], float)
    ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); vth = (x * A[:, :, 5] - y * A[:, :, 4]) / np.maximum(r, 0.5)
    C = vth[:, r < R - 1].mean(1); VL = vth[:, (R - r > 1) & (R - r < 2.5)].mean(1); tc = ts - (ts[1] - ts[0]) / 2
    t, fr = read_with_angle(d, np.array(Q[:2], float), 0.0); t = np.array(t); TT, ZZ, TM = [], [], []
    for tt, f in zip(t, fr):
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); w = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if w.sum() < N - 2: continue
        TT.append(np.degrees(np.angle(np.mean(np.exp(1j * (p[w, 3] - ph[w])))))); ZZ.append(np.mean(np.exp(1j * N * ph[w]))); TM.append(tt)
    TT, ZZ, TM = np.array(TT), np.array(ZZ), np.array(TM); X = np.unwrap(np.angle(ZZ)) / N * Rl
    def plat(arr, tt, a, b): m = (tt >= a) & (tt < b); return np.mean(arr[m]) if m.any() else np.nan
    C0, C1 = plat(C, tc, 1000, 1500), plat(C, tc, 2500, 3000); V0, V1 = plat(VL, tc, 1000, 1500), plat(VL, tc, 2500, 3000)
    T0, T1 = plat(TT, TM, 1000, 1500), plat(TT, TM, 2500, 3000); X0, X1 = plat(X, TM, 1000, 1500), plat(X, TM, 2500, 3000); dX = X1 - X0
    def zc(arr, tt, start):
        m = tt > start; a = arr[m]; b = tt[m]; s0 = np.sign(np.mean(arr[(tt > start - 300) & (tt < start)]))
        i = np.argmax(np.sign(a) == -s0) if (np.sign(a) == -s0).any() else None
        return b[i] - start if i is not None else np.nan
    tC = zc(C, tc, TF); tT = zc(np.convolve(TT, np.ones(5) / 5, "same"), TM, TF)
    win = (TM > TF - 100) & (TM < TF + 1400); Xs = np.convolve(X, np.ones(5) / 5, "same"); f = (Xs[win] - X0) / dX if abs(dX) > 0.3 else None
    if f is not None:
        tw = TM[win]; t10 = tw[np.argmax(f > 0.1)]; t90 = tw[np.argmax(f > 0.9)]; dur = t90 - t10
        vlay = abs(dX) / dur if dur > 0 else np.nan; vfl = np.mean(np.abs(VL[(tc >= t10) & (tc <= t90)])) if dur > 0 else np.nan
    else: t10 = t90 = dur = vlay = vfl = np.nan
    say(f"{name:22s}: circulation {C0:+.3f} -> {C1:+.3f} (zero crossing {tC:.0f} t.u. after the flip); fluid at layer {V0:+.3f} -> {V1:+.3f}; "
        f"tilt {T0:+.0f} -> {T1:+.0f} deg (zero crossing {tT:.0f} t.u.); registry shift {dX:+.2f} cells (lock: {np.radians(T1 - T0) * Rl / N:+.2f}), "
        f"transit {dur:.0f} t.u. from {t10 - TF:+.0f}; layer speed {vlay:.4f} vs fluid at layer {vfl:.3f} -> ratio {vlay / vfl if vfl else np.nan:.3f}")
    ax[row, 0].plot(tc, C, lw=0.7, label=name); ax[row, 1].plot(TM, TT, lw=0.5); ax[row, 2].plot(TM, X - X0, lw=0.6)
for i, lab in enumerate(("active", "passive")):
    ax[i, 0].set_ylabel(f"{lab}: circulation"); ax[i, 1].set_ylabel("wall tilt (deg)"); ax[i, 2].set_ylabel("registry shift (cells)")
    for a in ax[i]: a.axvline(TF, color="k", lw=0.6, ls="--"); a.set_xlabel("production time")
    ax[i, 0].legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(H, "swirl_response.png"), dpi=150); open(os.path.join(H, "swirl_response.txt"), "w").write("\n".join(L) + "\n")
