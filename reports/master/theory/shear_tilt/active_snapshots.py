#!/usr/bin/env python3
"""Snapshots of the active continuum disk (Q-tensor + backflow + activity, one constant; as active_spacing.py) for a movie (user, 2026-10-16: crowded vs
sparse layer, MPCD loses tilt, continuum does not). R 31.83, lambda 2, zeta 0.3, k 4 and 16 (2k-2 ring + bulk pair seeded); after 3000 t.u. the fields
are saved every 2 t.u. for 300 t.u. -> snaps_k{k}.npz (q1, q2, vx, vy per frame, times)."""
import sys, numpy as np
import shear_tilt as M
k = int(sys.argv[1]); M.k = k; M.qb1 = 0.5 * M.S0 * np.cos(2 * k * M.PH); M.qb2 = 0.5 * M.S0 * np.sin(2 * k * M.PH)
import active_bf as B
obs, (q1, q2, vx, vy) = B.run(2.0, 0.3, 0.0, ring=2 * k - 2, T0=300, T1=3000, dt=0.02, every=250, nflow=1)
lam_c = 2.0 * M.S0; dt = 0.02; Q1, Q2, VX, VY, T = [], [], [], [], []; i = 0; t = 0.0
fl = B.flowfields(vx, vy)
for n in range(150):
    for _ in range(int(2.0 / dt)):
        vx, vy = B.stokes(*B.forces(q1, q2, 0.3, lam_c, 0.0)); fl = B.flowfields(vx, vy); q1, q2 = M.step(q1, q2, dt, lam_c, fl)
    t += 2.0; Q1.append(q1.astype(np.float32)); Q2.append(q2.astype(np.float32)); VX.append(vx.astype(np.float32)); VY.append(vy.astype(np.float32)); T.append(t)
np.savez_compressed(f"snaps_k{k:02d}.npz", q1=np.array(Q1), q2=np.array(Q2), vx=np.array(VX), vy=np.array(VY), t=np.array(T), X=M.X, Y=M.Y, R=M.R)
print("done", k, obs[-1].get("tilt"), obs[-1].get("nwall"))
