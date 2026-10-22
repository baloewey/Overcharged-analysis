#!/usr/bin/env python3
"""Tracker biases on seeded layers with known truth (user, 2026-10-22). runs/rest_pinning: passive circle R 31.83, ring of 2k-2 wall +1/2 seeded
with the exact one-constant director (seed_director) at registry offset 0, 1/6, 1/3 of the spacing, i.e. true tilt tau0 = 0, 60, 120 deg and
true arg Z_N = tau0 (D_true = 0), k 4, 6, 10, 16, plus a bulk pair. The solver's tracker output in the first frames (t = 0 and 2 t.u., before
relaxation; the director is drawn from the seeded field with MPCD noise) vs the seeded positions (dirinit_seed.txt):
per wall defect, tangential and radial position error (cells), tilt error (tracker tilt - tau0); per frame, error of arg Z_N (from positions)
and of tau (layer mean), and of D = arg Z_N - tau. -> seeded_bias.txt"""
import glob, json, os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle
RD = "/home/baloewe/Physics/Overcharged/runs/rest_pinning"; L = ["run | spacing | frame t | matched | tangential err mean (sd) | radial err mean | tilt err mean (sd) deg | argZ err deg | tau err deg | D err deg"]
S = []
for run in sorted(glob.glob(RD + "/k*_t*")):
    nm = os.path.basename(run); k = int(nm[1:3]); tau0 = np.radians(int(nm[5:8])); N = 2 * k - 2
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; c = np.array(bc["Q"][:2])
    seed = np.loadtxt(run + "/dirinit_seed.txt"); sw = seed[(np.hypot(seed[:, 0], seed[:, 1]) > R - 6) & (seed[:, 2] > 0)]; phs = np.arctan2(sw[:, 1], sw[:, 0]); rs = np.hypot(sw[:, 0], sw[:, 1])
    zt = np.angle(np.exp(1j * N * phs).mean())
    t, fr = read_with_angle(run, c, 0.0)
    for ti, f in zip(t, fr):
        if ti > 2.01: break
        q = f[f[:, 2] > 0]; ph = np.arctan2(q[:, 1], q[:, 0]); r = np.hypot(q[:, 0], q[:, 1]); m = R - r < 6; q, ph, r = q[m], ph[m], r[m]
        dt, dr, dtil = [], [], []
        for j in range(len(sw)):
            dphi = np.angle(np.exp(1j * (ph - phs[j]))); i = np.argmin(np.abs(dphi) * rs[j] + np.abs(r - rs[j]))
            if abs(dphi[i]) * rs[j] > 3: continue
            dt.append(dphi[i] * rs[j]); dr.append(r[i] - rs[j]); dtil.append(np.angle(np.exp(1j * (q[i, 3] - ph[i] - tau0))))
        if len(q) == N:
            za = np.angle(np.exp(1j * N * ph).mean()); tm = np.angle(np.mean(np.exp(1j * (q[:, 3] - ph))))
            ez, et = np.angle(np.exp(1j * (za - zt))), np.angle(np.exp(1j * (tm - tau0))); ed = np.angle(np.exp(1j * (ez - et)))
        else: ez = et = ed = np.nan
        dt, dr, dtil = map(np.array, (dt, dr, dtil)); sp = 2 * np.pi * rs.mean() / N
        S.append((k, sp, int(nm[5:8]), ti, np.mean(dt), np.mean(dr), np.degrees(np.mean(dtil)), np.degrees(ez), np.degrees(et), np.degrees(ed)))
        L.append(f"{nm} | {sp:5.1f} | {ti:3.0f} | {len(dt)}/{N} | {np.mean(dt):+.3f} ({np.std(dt):.2f}) | {np.mean(dr):+.3f} | {np.degrees(np.mean(dtil)):+6.1f} ({np.degrees(np.std(dtil)):5.1f}) | " +
                 (f"{np.degrees(ez):+6.1f} | {np.degrees(et):+6.1f} | {np.degrees(ed):+6.1f}" if np.isfinite(ez) else f"excluded: tracker found {len(q)} wall +1/2, not {N}"))
S = np.array(S); L.append("\nBY k (mean over tau0 and the two frames): argZ err, tau err, D err (deg); tangential position err (cells)")
for k in (4, 6, 10, 16):
    m = (S[:, 0] == k) & np.isfinite(S[:, 9]); 
    if m.any(): L.append(f"k {k:2d} (spacing {S[m,1].mean():.1f}): argZ {np.mean(S[m,7]):+5.1f} +- {np.std(S[m,7]):.1f}, tau {np.mean(S[m,8]):+5.1f} +- {np.std(S[m,8]):.1f}, D {np.mean(S[m,9]):+5.1f} +- {np.std(S[m,9]):.1f}; tangential {np.mean(S[m,4]):+.3f} cells")
open("seeded_bias.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
