#!/usr/bin/env python3
"""Is there a flow-alignment signature near the wall? Measured director minus the one-constant director of the measured defects.
Elastic field: theta_el(z) = sum_i q_i [arg(z - z_i) + arg(1 - z conj(z_i)/R^2)] (all +-1/2 of the frame; anchoring theta = k phi on the wall,
exact). Residual delta = theta_meas - theta_el (mod pi, in (-pi/2, pi/2]) on the cell grid, signed by the circulation of the flow interval ending
at that time (s = sign <v_theta>), averaged in bands of distance from the wall; also correlated with the local shear rate
d v_theta / dr computed from the interval-averaged flow. Flow alignment (lambda = 2, Leslie angle 30 deg from the flow) predicts a residual
of fixed sign relative to the shear, growing with |shear|, largest where the elastic torque is weakest.
Runs: active circles k = 2, 6, 10, R = 31.83 (Rscan_R31.83_seeds/s02), director snapshots every 50 t.u. with matching defect frames.
Writes residual_director.txt."""
import gzip, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
def frames(fn, ncol):
    rows = np.loadtxt((l for l in gzip.open(fn, "rt") if len(l.split()) == ncol and l.split()[0][0].isdigit())); return rows
def defects(fn):
    F = {}; cur = None
    for l in gzip.open(fn, "rt"):
        p = l.split()
        if len(p) == 2:
            try: cur = round(float(p[0]), 3); F[cur] = []
            except ValueError: cur = None
        elif len(p) == 4 and cur is not None:
            try: F[cur].append([float(x) for x in p])
            except ValueError: pass
    return F
for k in (2, 6, 10):
    run = f"{ROOT}/cavity/Rscan_R31.83_seeds/s02/k{k:02d}"; bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; R = bc["R"]
    Dr = frames(run + "/directorfield.dat.gz", 8); Fl = frames(run + "/flowfield.dat.gz", 7); Df = defects(run + "/defects.dat.gz")
    td = np.unique(Dr[:, 0]); td = td[td > td.max() * 0.2]
    bands = np.arange(0, 12, 1.0); acc = {b: [] for b in bands}; shear_res = []
    for t in td:
        if round(t, 3) not in Df or not Df[round(t, 3)]: continue
        d = Dr[Dr[:, 0] == t]; f = Fl[np.isclose(Fl[:, 0], t)]
        if len(f) == 0: continue
        x = d[:, 1] + 0.5 - Q[0]; y = d[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); m = (r < R - 0.3) & (d[:, 7] > 0.3)
        th_m = np.arctan2(d[m, 5], d[m, 4]); z = x[m] + 1j * y[m]
        th_el = np.zeros(m.sum())
        for xd, yd, q, ang in Df[round(t, 3)]:
            zi = (xd - Q[0]) + 1j * (yd - Q[1])
            th_el += q * (np.angle(z - zi + 1e-9) + np.angle(1 - z * np.conj(zi) / R ** 2))
        delta = np.angle(np.exp(2j * (th_m - th_el))) / 2
        fx = f[:, 1] + 0.5 - Q[0]; fy = f[:, 2] + 0.5 - Q[1]; fr = np.hypot(fx, fy); vt = (fx * f[:, 5] - fy * f[:, 4]) / np.maximum(fr, 0.5)
        ins = fr < R - 1; s = np.sign(np.mean(vt[ins]))
        if s == 0: continue
        dep = R - r[m]
        for b in bands:
            sel = (dep >= b) & (dep < b + 1)
            if sel.any(): acc[b].append(s * np.mean(delta[sel]))
        # local shear rate: radial profile of v_theta (2-cell bins), derivative at each director cell's radius
        rb = np.arange(0, R, 1.0); prof = np.array([np.mean(vt[(fr >= a) & (fr < a + 1)]) if ((fr >= a) & (fr < a + 1)).any() else np.nan for a in rb])
        dv = np.gradient(prof, 1.0) - prof / np.maximum(rb + 0.5, 0.5)            # r d(v/r)/dr: the shear rate of a swirling flow
        sh = np.interp(r[m], rb + 0.5, np.nan_to_num(dv))
        near = dep < 6
        shear_res.append((sh[near], delta[near]))
    say(f"k={k}: {len(acc[0])} frames; signed mean residual (deg, sign of the circulation) by distance from the wall: " +
        " ".join(f"{b:.0f}-{b + 1:.0f}:{np.degrees(np.mean(acc[b])):+.1f}" for b in bands if acc[b]))
    shv = np.concatenate([a for a, _ in shear_res]); dlv = np.concatenate([b for _, b in shear_res])
    q = np.quantile(shv, [0, 0.2, 0.4, 0.6, 0.8, 1])
    say("      residual vs local shear r d(v/r)/dr (within 6 cells of the wall): " + " ".join(
        f"[{q[i]:+.3f},{q[i + 1]:+.3f}]:{np.degrees(np.mean(dlv[(shv >= q[i]) & (shv <= q[i + 1])])):+.1f}" for i in range(5)) +
        f"; corr {np.corrcoef(shv, dlv)[0, 1]:+.2f}")
open(os.path.join(H, "residual_director.txt"), "w").write("\n".join(L) + "\n")
