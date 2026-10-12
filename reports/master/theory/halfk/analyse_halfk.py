#!/usr/bin/env python3
"""Half-integer windings (runs/halfk: k = 0.5, 1.5, 2.5, 3.5, 5.5; R = 31.83; active act 0.07 x2 seeds, passive x1; production binary, eMode 0).
The wall director theta = k phi is (k - 1) phi from the normal: homeotropic points phi_n = 2 pi n/(2k - 2) (an odd number, 2k - 2, for half-integer k;
none for k = 1, and for k = 0.5 one 'anti' point per turn). Per run, second half:
  counts: wall (depth < 4) and bulk (depth > 6) +1/2 and -1/2 (means, modal wall +1/2 count); bulk charge Q_b; expected wall +1/2: 2k (passive), 2k - 2 (active)
  bulk +1: charge-weighted centroid of the bulk defects, distance from the centre (/R); for k = 1.5 the angle between it and the single wall +1/2
  registry: frames with the modal wall count N_w: Z = <exp(i N_w phi)>, |Z| and arg Z; lock D = arg Z - tau (tau = mean tail tilt) and its spread
  tilt |tau| and coherence; circulation (flow output, active) |C| and sign changes
Reference: integer k = 2..6 at the same R (Rscan_R31.83_seeds/s02). Writes analyse_halfk.txt."""
import glob, gzip, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
def analyse(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; R = bc["R"]; Q = bc["Q"]; k = bc["customAnch"]
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5)
    wp, wm, bp, bm, qb, cen, ang15, Zs, taus, nws = [], [], [], [], [], [], [], [], [], []
    for f in fr:
        r = np.hypot(f[:, 0], f[:, 1]); dep = R - r; ph = np.arctan2(f[:, 1], f[:, 0]); pl = f[:, 2] > 0
        w = dep < 4; b = dep > 6
        wp.append((w & pl).sum()); wm.append((w & ~pl).sum()); bp.append((b & pl).sum()); bm.append((b & ~pl).sum())
        q = 0.5 * np.sign(f[b, 2]).sum(); qb.append(q)
        if b.any() and abs(q) > 0.1:
            c = np.sum(0.5 * np.sign(f[b, 2]) * (f[b, 0] + 1j * f[b, 1])) / q; cen.append(c)
            if abs(k - 1.5) < 1e-9 and (w & pl).sum() == 1: ang15.append(np.degrees(np.angle(c * np.exp(-1j * ph[w & pl][0]))))
        nws.append((w & pl).sum())
        if (w & pl).any():
            Zs.append(((w & pl).sum(), np.mean(np.exp(1j * (w & pl).sum() * ph[w & pl])), np.angle(np.mean(np.exp(1j * (f[w & pl, 3] - ph[w & pl]))))))
    nws = np.array(nws); mode = np.bincount(nws).argmax()
    sel = [(z, ta) for n, z, ta in Zs if n == mode]
    Z = np.mean([z for z, _ in sel]) if sel else np.nan; tau = np.array([ta for _, ta in sel])
    D = np.array([np.angle(z * np.exp(-1j * ta)) for z, ta in sel]) if sel else np.array([np.nan])
    out = dict(k=k, R=R, wp=np.mean(wp), wm=np.mean(wm), bp=np.mean(bp), bm=np.mean(bm), qb=np.mean(qb), mode=mode, fmode=np.mean(nws == mode),
               cen=np.mean(np.abs(cen)) / R if cen else np.nan, ang15=(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(ang15))))), abs(np.mean(np.exp(1j * np.radians(ang15))))) if ang15 else None,
               Zabs=abs(Z) if sel else np.nan, argZ=np.angle(Z) if sel else np.nan, Dmean=np.angle(np.mean(np.exp(1j * D))), Dconc=abs(np.mean(np.exp(1j * D))),
               tilt=np.degrees(np.mean(np.abs(tau))) if len(tau) else np.nan, coh=abs(np.mean(np.exp(1j * tau))) if len(tau) else np.nan)
    ff = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else (run + "/flowfield.dat" if os.path.exists(run + "/flowfield.dat") else None)
    if ff:
        op = gzip.open if ff.endswith(".gz") else open
        A = np.array([l.split() for l in op(ff, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()], float)
        ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7); A = A[len(ts) // 2:]
        x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); m = r < R - 1
        C = ((x * A[:, :, 5] - y * A[:, :, 4]) / np.maximum(r, 0.5))[:, m].mean(1)
        out.update(circ=np.mean(np.abs(C)), flips=int(np.sum(np.diff(np.sign(C)) != 0)), nint=len(C))
    return out
hdr = f"{'run':24s} {'k':>4s} | {'wall +':>6s} {'wall -':>6s} {'bulk +':>6s} {'bulk -':>6s} {'Q_bulk':>6s} | {'mode':>4s} {'P':>4s} (exp.) | {'bulk centroid r/R':>8s} | {'|Z|':>5s} {'argZ':>6s} {'D':>6s} {'D conc':>6s} | {'|tilt|':>6s} {'coh':>5s} | {'|circ|':>6s} {'flips':>5s}"
for title, runs in (("half-integer, active", sorted(glob.glob(ROOT + "/halfk/active_*"))), ("half-integer, passive", sorted(glob.glob(ROOT + "/halfk/passive_*"))),
                    ("integer reference (active, R 31.83, s02)", [ROOT + f"/cavity/Rscan_R31.83_seeds/s02/k{k:02d}" for k in (2, 3, 4, 5, 6)])):
    say(f"\n== {title}"); say(hdr)
    for run in runs:
        if not os.path.exists(run + "/DONE"): say(f"{run}: not finished"); continue
        o = analyse(run); passive = "passive" in run; exp_ = 2 * o["k"] if passive else 2 * o["k"] - 2
        say(f"{run.split('runs/')[1][-24:]:24s} {o['k']:4.1f} | {o['wp']:6.2f} {o['wm']:6.2f} {o['bp']:6.2f} {o['bm']:6.2f} {o['qb']:+6.2f} | {o['mode']:4d} {o['fmode']:4.2f} ({exp_:+.0f}) | "
            f"{o['cen']:8.2f} | {o['Zabs']:5.2f} {o['argZ']:+6.2f} {o['Dmean']:+6.2f} {o['Dconc']:6.2f} | {o['tilt']:6.1f} {o['coh']:5.2f} | {o.get('circ', np.nan):6.3f} {o.get('flips', -1):5d}"
            + (f" | k=1.5: bulk centroid at {o['ang15'][0]:+.0f} deg from the wall defect (concentration {o['ang15'][1]:.2f})" if o['ang15'] else ""))
open(os.path.join(H, "analyse_halfk.txt"), "w").write("\n".join(L) + "\n")
