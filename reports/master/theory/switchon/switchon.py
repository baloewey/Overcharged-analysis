#!/usr/bin/env python3
"""Passive -> active switch (runs/switchon; activity 0 -> 0.07 at t = 3000, k = 6, R = 31.83, 8 initial conditions; 2026-10-22).
Per run, in 20 t.u. bins: wall +1/2 count (within 4 cells of the wall), tilt tau = circular mean of (tail angle - phi) over wall +1/2, coherence |<e^{i tau}>|,
registry arg Z_10 (Z_10 = <exp(10 i phi)> over wall +1/2) and lock-corrected D = arg Z_10 - tau (frames with exactly 10 wall defects),
arg Z_12 (passive ring), and circulation <v_phi>/<|v|> from the flow field (every 10 t.u.). Passive registry arg Z_12 from the second half of the passive
runs/pinning_ic run (frames with exactly 12 wall +1/2).
Final state (last 1000 t.u.) compared with runs/pinning_ic_active (active from the start: arg Z_10 = +-1.17, D = 0.03 +- 0.07).
-> switchon.txt, switchon.npz"""
import glob, gzip, json, os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle
RD = "/home/baloewe/Physics/Overcharged/runs/switchon"; BIN = 20.0; out, S = [], {}
def flow_circ(run, c, R):
    T, C = [], []; cur_t, num, den = None, 0.0, 0.0
    with gzip.open(run + "/flowfield.dat.gz", "rt") as f:
        for line in f:
            p = line.split()
            if len(p) != 7 or not p[0][0].isdigit(): continue
            t = float(p[0])
            if cur_t is not None and t != cur_t: T.append(cur_t); C.append(num / max(den, 1e-12)); num = den = 0.0
            cur_t = t; x, y = float(p[1]) + 0.5 - c[0], float(p[2]) + 0.5 - c[1]; r = np.hypot(x, y)
            if r < 0.5 or r > R - 1: continue
            vx, vy = float(p[4]), float(p[5]); num += (-vx * y + vy * x) / r; den += np.hypot(vx, vy)
    if cur_t is not None: T.append(cur_t); C.append(num / max(den, 1e-12))
    return np.array(T), np.array(C)
for run in sorted(glob.glob(RD + "/*/")):
    run = run.rstrip("/"); nm = os.path.basename(run); bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; c = np.array(bc["Q"][:2])
    t, fr = read_with_angle(run, c, 0.0); t = np.array(t)
    nw, tau, coh, z10, d10, z12 = [], [], [], [], [], []
    for f in fr:
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); m = R - np.hypot(p[:, 0], p[:, 1]) < 4; nw.append(m.sum())
        if m.sum() == 0: tau.append(np.nan); coh.append(np.nan); z10.append(np.nan); d10.append(np.nan); z12.append(np.nan); continue
        e = np.exp(1j * (p[m, 3] - ph[m])); tau.append(np.angle(e.mean())); coh.append(abs(e.mean()))
        z12.append(np.angle(np.exp(12j * ph[m]).mean()) if m.sum() == 12 else np.nan)
        if m.sum() == 10: z = np.exp(10j * ph[m]).mean(); z10.append(np.angle(z)); d10.append(np.angle(np.exp(1j * (np.angle(z) - tau[-1]))))
        else: z10.append(np.nan); d10.append(np.nan)
    nw, tau, coh, z10, d10, z12 = map(np.array, (nw, tau, coh, z10, d10, z12)); tf, cf = flow_circ(run, c, R)
    S[nm] = dict(t=t, nw=nw, tau=tau, coh=coh, z10=z10, d10=d10, z12=z12, tf=tf, circ=cf)
    # event times
    def first(cond):
        i = np.flatnonzero(cond); return t[i[0]] if len(i) else np.nan
    t10 = first((t > 3000) & (nw <= 10)); tcoh = first((t > 3000) & (coh > 0.9)); tc = tf[np.flatnonzero((tf > 3000) & (np.abs(cf) > 0.3))[0]] if np.any((tf > 3000) & (np.abs(cf) > 0.3)) else np.nan
    late = t > t.max() - 1000; m10 = late & np.isfinite(d10)
    tp, fp = read_with_angle("/home/baloewe/Physics/Overcharged/runs/pinning_ic/" + nm, c, 0.5); zp = []   # passive registry: second half of the passive run
    for f in fp:
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); m = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if m.sum() == 12: zp.append(np.exp(12j * ph[m]).mean())
    pz12 = np.angle(np.mean(np.array(zp) / np.abs(zp))) if zp else np.nan
    Dm = np.mean(np.exp(1j * d10[m10])) if m10.any() else np.nan; Am = np.mean(np.exp(1j * z10[m10])) if m10.any() else np.nan
    lc = late & (np.abs(tf[:0]).size == 0)
    cl = np.nanmean(cf[tf > tf.max() - 1000])
    out.append(f"{nm:9s} | passive arg Z12 {pz12:+.2f} ({len(zp)} frames) | first 10 wall {t10 - 3000:6.0f} t.u. after switch | coherent (C>0.9) {tcoh - 3000:6.0f} | |circ|>0.3 {tc - 3000:6.0f} | "
               f"last 1000 t.u.: wall {np.mean(nw[late]):.2f}, frames with 10: {m10.mean() / late.mean():.2f}, tilt {np.degrees(np.angle(np.mean(np.exp(1j * tau[late] * np.sign(np.nanmean(tau[late])))))) * np.sign(np.nanmean(tau[late])):+.1f}, "
               f"arg Z10 {np.angle(Am):+.2f} (|.| {abs(Am):.2f}), D {np.angle(Dm):+.2f} (|.| {abs(Dm):.2f}), circ {cl:+.2f}")
np.savez("switchon.npz", **{f"{k}_{q}": v for k, d in S.items() for q, v in d.items()})
open("switchon.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
