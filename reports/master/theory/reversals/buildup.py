#!/usr/bin/env python3
"""Tilt build-up from a random start (runs/tilt_buildup, k = 6, R = 31.83, act 0.07, no warm-up): how does the wall layer respond to the flow?
Per run: wall layer (depth < 4, >= 2k - 4 wall +1/2): tilt tau(t) (tail convention, signed by its final sign), registry shift X(t) =
unwrap(arg <exp(i N phi)>) R_l / N (cells, signed likewise), fluid tangential velocity at the layer depth (1-2.5 cells) per 10 t.u.
Reported: tau and X at t = 0, 25, 50, 100, 200, 400, 600; the time to reach half and 90% of the final tilt; the lock (slope of X vs the predicted
tau R_l/(2k - 2)); the effective advection coefficient beta = <dX/dt> / <v_fluid at layer> over the build-up (until 90% of the final tilt) and
in the steady state. Writes buildup.txt and buildup.png."""
import glob, gzip, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for d in sorted(glob.glob(ROOT + "/tilt_buildup/s*")):
    if not os.path.exists(d + "/DONE"): say(f"{d}: not finished"); continue
    bc = json.load(open(d + "/input.json"))["BC"][0]; R = bc["R"]; Q = bc["Q"]; k = int(bc["customAnch"]); N = 2 * k - 2; Rl = R - 1.6
    t, fr = read_with_angle(d, np.array(Q[:2], float), 0.0); t = np.array(t)
    fn = d + "/flowfield.dat.gz" if os.path.exists(d + "/flowfield.dat.gz") else d + "/flowfield.dat"; op = gzip.open if fn.endswith(".gz") else open
    A = np.array([l.split() for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()], float)
    ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); lay = (R - r > 1.0) & (R - r < 2.5)
    VL = ((x * A[:, :, 5] - y * A[:, :, 4]) / np.maximum(r, 0.5))[:, lay].mean(1); tv = ts - (ts[1] - ts[0]) / 2
    TT, ZZ, TM = [], [], []
    for tt, f in zip(t, fr):
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); w = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if w.sum() < N - 2: continue
        TT.append(np.angle(np.mean(np.exp(1j * (p[w, 3] - ph[w]))))); ZZ.append(np.mean(np.exp(1j * N * ph[w]))); TM.append(tt)
    TT, ZZ, TM = np.array(TT), np.array(ZZ), np.array(TM)
    if len(TM) < 20: say(f"{d}: too few layer frames ({len(TM)})"); continue
    X = np.unwrap(np.angle(ZZ)) / N * Rl
    late = TM > TM.max() * 0.6; sg = np.sign(np.mean(np.sin(TT[late]))) or 1.0
    tau = np.degrees(TT) * sg; X = (X - X[0]) * sg; Vs = VL * sg
    tf = np.mean(tau[late]); sm = np.convolve(tau, np.ones(5) / 5, "same")
    t50 = TM[np.argmax(sm > 0.5 * tf)] if (sm > 0.5 * tf).any() else np.nan; t90 = TM[np.argmax(sm > 0.9 * tf)] if (sm > 0.9 * tf).any() else np.nan
    pred = np.radians(tau) * Rl / N; Xc = X - X[late].mean() + pred[late].mean()
    slope = np.polyfit(pred, Xc, 1)[0]
    bu = TM <= t90 if np.isfinite(t90) else TM <= TM.max()
    dX = (X[bu][-1] - X[bu][0]) / max(TM[bu][-1] - TM[bu][0], 1); vfl = np.mean(Vs[tv <= (t90 if np.isfinite(t90) else tv.max())])
    dXs = (X[late][-1] - X[late][0]) / max(TM[late][-1] - TM[late][0], 1); vfs = np.mean(Vs[tv > TM.max() * 0.6])
    say(f"{os.path.basename(d)}: first layer frame t={TM[0]:.0f}; tilt (deg) at t=" + " ".join(f"{q}:{tau[np.argmin(np.abs(TM - q))]:+.0f}" for q in (0, 25, 50, 100, 200, 400, 600)) +
        f"; final {tf:+.0f}; t_50 {t50:.0f}, t_90 {t90:.0f} t.u.")
    say(f"     shift X (cells) at t=" + " ".join(f"{q}:{X[np.argmin(np.abs(TM - q))]:+.2f}" for q in (0, 25, 50, 100, 200, 400, 600)) +
        f"; lock slope X vs tau R_l/N {slope:.2f}; build-up: layer {dX:+.4f} vs fluid at layer {vfl:+.3f} cells/t.u. -> beta {dX / vfl if vfl else np.nan:+.3f}; "
        f"steady: layer {dXs:+.4f} vs fluid {vfs:+.3f} -> beta {dXs / vfs if vfs else np.nan:+.3f}")
    ax[0].plot(TM, tau, lw=0.7); ax[1].plot(TM, X, lw=0.7); ax[2].plot(tv, Vs, lw=0.7)
ax[0].set_ylabel("wall tilt (deg, signed)"); ax[1].set_ylabel("registry shift X (cells)"); ax[2].set_ylabel("fluid speed at layer depth")
for a in ax: a.set_xlabel("t (from a random start)")
fig.tight_layout(); fig.savefig(os.path.join(H, "buildup.png"), dpi=150); open(os.path.join(H, "buildup.txt"), "w").write("\n".join(L) + "\n")
