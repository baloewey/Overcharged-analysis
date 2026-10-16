#!/usr/bin/env python3
"""Is the tilt measured the same way everywhere? (user, 2026-10-16) Check of the MPCD tracker angle (defects.dat 'angle', used by every MPCD tilt
analysis) against the direction of div Q (the comet's radial arm; used by the continuum codes) on SEEDED fields with known geometry
(runs/rest_pinning: passive, seeded 2k-2 ring at registry offsets 0, 1/6, 1/3 -> intended tilt 0, 60, 120 deg by the lock). At t = 0 the particle
orientations are the seeded field, so tracker angle and div Q (from dirinit.dat, averaged over a disc of radius r_d around the seed) can be compared
defect by defect. -> tracker_check.txt"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); RP = "/home/baloewe/Physics/Overcharged/runs/rest_pinning"
L = ["TRACKER ANGLE vs div Q on seeded rings (t = 0). Tilt = angle from the outward normal (deg), median over wall defects; spread = sd of tracker - divQ.",
     "run | spacing | intended tilt (lock) | tilt from div Q (disc r 1.5 / 2 / 3) | tracker tilt | tracker - divQ(r 2): median, spread | seed-tracker distance"]
for run in sorted(glob.glob(RP + "/k*_t*")):
    bc = json.load(open(run + "/input.json"))["BC"][0]; Q = np.array(bc["Q"][:2]); R = bc["R"]; k = bc["customAnch"]; tag = int(run[-3:])
    if not os.path.exists(run + "/defects.dat") and not os.path.exists(run + "/defects.dat.gz"): continue
    D = np.loadtxt(run + "/dirinit.dat"); n = int(D[:, 0].max()) + 1; m_ = int(D[:, 1].max()) + 1; th = np.full((n, m_), np.nan); th[D[:, 0].astype(int), D[:, 1].astype(int)] = D[:, 2]
    q1, q2 = np.nan_to_num(np.cos(2 * th)), np.nan_to_num(np.sin(2 * th))
    dvx = np.gradient(q1, axis=0) + np.gradient(q2, axis=1); dvy = np.gradient(q2, axis=0) - np.gradient(q1, axis=1)
    X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m_) + 0.5 - Q[1], indexing="ij")
    seeds = np.loadtxt(run + "/dirinit_seed.txt"); t, fr = read_with_angle(run, Q, 0.0); f0 = fr[0]; rows = []
    for xs, ys, qq in seeds:
        if qq < 0 or np.hypot(xs, ys) < R - 6: continue
        phi = np.arctan2(ys, xs); w = lambda a: np.degrees(np.angle(np.exp(1j * a)))
        dq = [w(np.arctan2(dvy[np.hypot(X - xs, Y - ys) < rd].mean(), dvx[np.hypot(X - xs, Y - ys) < rd].mean()) - phi) for rd in (1.5, 2.0, 3.0)]
        j = np.argmin(np.hypot(f0[:, 0] - xs, f0[:, 1] - ys)); rows.append(dq + [w(f0[j, 3] - phi), w(np.radians(w(f0[j, 3] - phi) - dq[1])), np.hypot(f0[j, 0] - xs, f0[j, 1] - ys)])
    o = np.array(rows)
    L.append(f"{os.path.basename(run)} | {2*np.pi*R/(2*k-2):5.1f} | {tag:4d} | {np.median(o[:,0]):6.1f} {np.median(o[:,1]):6.1f} {np.median(o[:,2]):6.1f} | {np.median(o[:,3]):6.1f} | {np.median(o[:,4]):+6.1f} {np.std(o[:,4]):5.1f} | {np.median(o[:,5]):.2f}")
print("\n".join(L)); open(os.path.join(H, "tracker_check.txt"), "w").write("\n".join(L) + "\n")
