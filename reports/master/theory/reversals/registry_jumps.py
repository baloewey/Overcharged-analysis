#!/usr/bin/env python3
"""Registry jumps of the wall layer at circulation reversals (activity 0.05 circles, k = 9, 10; H3pilot + seeds; flow every 50 t.u.).
Per run: circulation C(t) (interval means, assigned to t - 25), zero crossings = reversals. Wall layer (depth < 4, frames with at least 2k - 6 wall +1/2; the layer at activity 0.05 has vacancies):
  tilt tau(t) = circular mean of (tail angle - phi) [deg]; registry shift X(t) = unwrap(arg <exp(i N phi)>) R_l / N (cells; vacancy-tolerant); fluid tangential velocity at the layer depth (1-2.5 cells from the wall) per flow interval.
Per reversal (window +-400 t.u.): plateau values of tau and X before/after (means over -400..-150 and +150..+400), shift dX, predicted shift
2 tau R_l/(2k-2) (lock), jump duration (10-90% of dX), midpoint time relative to the circulation zero crossing, peak layer speed, mean fluid
speed at the layer depth during the jump, and their ratio (effective advection coefficient while the layer moves); sign of dX vs the new flow.
Writes registry_jumps.txt and registry_jumps.png."""
import glob, gzip, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
runs = sorted(set(glob.glob(ROOT + "/cavity/H3pilot/act0.05/k09") + glob.glob(ROOT + "/cavity/H3pilot/act0.05/k10") +
                  glob.glob(ROOT + "/cavity/H3pilot_seeds/act0.05/s0*/k09") + glob.glob(ROOT + "/cavity/H3pilot_seeds/act0.05/s0*/k10")))
events = []; fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for d in runs:
    bc = json.load(open(d + "/input.json"))["BC"][0]; R = bc["R"]; Q = bc["Q"]; k = int(bc["customAnch"]); N = 2 * k - 2
    op = gzip.open(d + "/flowfield.dat.gz", "rt") if os.path.exists(d + "/flowfield.dat.gz") else open(d + "/flowfield.dat")
    A = np.array([l.split() for l in op if len(l.split()) == 7 and l.split()[0][0].isdigit()], float)
    ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ok = (r < R - 1) & (r > 0.5); lay = (R - r > 1.0) & (R - r < 2.5)
    vth = (x * A[:, :, 5] - y * A[:, :, 4]) / np.maximum(r, 0.5)
    C = vth[:, ok].mean(1); VL = vth[:, lay].mean(1); tc = ts - 25
    t, fr = read_with_angle(d, np.array(Q[:2], float), 0.0); t = np.array(t)
    TT, ZZ, TM = [], [], []
    for tt, f in zip(t, fr):
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); w = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if w.sum() < N - 4: continue
        TT.append(np.degrees(np.angle(np.mean(np.exp(1j * (p[w, 3] - ph[w])))))); ZZ.append(np.mean(np.exp(1j * N * ph[w]))); TM.append(tt)
    TT, ZZ, TM = np.array(TT), np.array(ZZ), np.array(TM)
    if len(TM) < 100: continue
    Rl = R - 1.6; X = (np.unwrap(np.angle(ZZ)) - np.angle(ZZ[0])) / N * Rl          # registry shift along the wall (cells), vacancy-tolerant
    zc = [tc[i] + (tc[i + 1] - tc[i]) * C[i] / (C[i] - C[i + 1]) for i in range(len(C) - 1) if np.sign(C[i]) != np.sign(C[i + 1]) and np.sign(C[i]) != 0]
    for t0 in zc:
        pre = (TM > t0 - 400) & (TM < t0 - 150); post = (TM > t0 + 150) & (TM < t0 + 400); win = (TM > t0 - 400) & (TM < t0 + 400)
        if pre.sum() < 20 or post.sum() < 20: continue
        tau0, tau1 = TT[pre].mean(), TT[post].mean(); X0, X1 = X[pre].mean(), X[post].mean(); dX = X1 - X0
        Xw = X[win]; Tw = TM[win]; Xs = np.convolve(Xw, np.ones(5) / 5, "same")
        f = (Xs - X0) / dX if abs(dX) > 0.3 else np.full(len(Xs), np.nan)
        try:
            t10 = Tw[np.argmax(f > 0.1)]; t90 = Tw[np.argmax(f > 0.9)]; tmid = Tw[np.argmax(f > 0.5)]
        except Exception: t10 = t90 = tmid = np.nan
        dur = t90 - t10
        vlay = np.gradient(Xs, Tw); jm = (Tw >= t10) & (Tw <= t90)
        vpeak = np.max(np.abs(vlay[jm])) if jm.any() else np.nan; vmean = abs(dX) / dur if dur > 0 else np.nan
        fl = (tc >= t10 - 25) & (tc <= t90 + 25); vfl = np.mean(np.abs(VL[fl])) if fl.any() else np.nan
        newflow = np.sign(C[np.argmin(np.abs(tc - (t0 + 200)))])
        pred = (np.radians(tau1) - np.radians(tau0)) * Rl / N
        events.append(dict(run=d.split("cavity/")[1], k=k, t0=t0, tau0=tau0, tau1=tau1, dX=dX, pred=pred, dur=dur, tmid=tmid - t0, vmean=vmean, vpeak=vpeak, vfl=vfl,
                           ratio=vmean / vfl if vfl and vfl > 0 else np.nan, dir_vs_flow=np.sign(dX) * newflow))
        e = events[-1]
        say(f"{e['run']:30s} t0 {t0:7.0f}: tilt {tau0:+5.0f} -> {tau1:+5.0f} deg; shift {dX:+6.2f} cells (lock predicts {pred:+6.2f}); duration {dur:5.0f} t.u.; "
            f"midpoint {e['tmid']:+5.0f} t.u. after the flow reversal; layer speed {vmean:.4f} (peak {vpeak:.3f}) vs fluid at layer depth {vfl:.3f} "
            f"-> ratio {e['ratio']:.3f}; shift along the new flow: {'yes' if e['dir_vs_flow'] > 0 else 'no'}")
        sel = win
        ax[0].plot(TM[sel] - t0, TT[sel], lw=0.6); ax[1].plot(TM[sel] - t0, X[sel] - X0, lw=0.6)
if events:
    E = {kk: np.array([e[kk] for e in events], float) for kk in ("dX", "pred", "dur", "tmid", "ratio", "dir_vs_flow", "tau0", "tau1")}
    say(f"\n{len(events)} reversals: shift vs lock prediction: corr {np.corrcoef(E['dX'], E['pred'])[0, 1]:+.2f}, slope {np.polyfit(E['pred'], E['dX'], 1)[0]:.2f}; "
        f"|shift| median {np.median(np.abs(E['dX'])):.2f} cells (predicted {np.median(np.abs(E['pred'])):.2f}); duration median {np.nanmedian(E['dur']):.0f} t.u.; "
        f"midpoint median {np.nanmedian(E['tmid']):+.0f} t.u. after the flow reversal; advection ratio median {np.nanmedian(E['ratio']):.3f} "
        f"(10-90%: {np.nanpercentile(E['ratio'], 10):.3f}-{np.nanpercentile(E['ratio'], 90):.3f}); shift along the new flow in {np.mean(E['dir_vs_flow'] > 0):.0%}")
    ax[2].plot(E["pred"], E["dX"], "o"); ax[2].axline((0, 0), slope=1, color="k", lw=0.6)
ax[0].set_xlabel("t - t_reversal"); ax[0].set_ylabel("wall tilt (deg)"); ax[1].set_xlabel("t - t_reversal"); ax[1].set_ylabel("layer shift (cells)")
ax[2].set_xlabel("lock prediction (cells)"); ax[2].set_ylabel("measured shift (cells)")
for a, tl in zip(ax, ("(a) tilt across reversals", "(b) collective layer displacement", "(c) shift vs tilt-lock prediction")): a.set_title(tl, fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(H, "registry_jumps.png"), dpi=150)
open(os.path.join(H, "registry_jumps.txt"), "w").write("\n".join(L) + "\n")
