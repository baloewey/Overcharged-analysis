#!/usr/bin/env python3
"""Stiffness of the flow-selected registry (user list, 2026-10-22). In the lock + flow picture the collective registry shift D (all wall defects
moved along the wall) rotates the far field by kappa D (kappa = pi/a, a = site spacing); the shear boundary layer resists a far-field rotation
delta with a torque per unit length K delta / l_E, l_E = sqrt(K / (chi gamma1 s sqrt(lambda^2-1))). The generalised force on D per defect is
then K kappa^2 a D / l_E, i.e. a spring k_D = K pi^2 / (a l_E) per defect, N k_D for the common shift; equipartition at an effective
temperature T: rms(D_c) = sqrt(T / (N k_D)).
Measured: passive circles under an imposed vortex (runs/swirl_crowd, swirl_shear; thermal noise only, T_eff should be ~1) and active circles of the
main sweep (runs/cavity/server_main/main_act0.07, k 3-10, seeds 1-3); spacing a from tiltcrowd/defect_flow, wall shear s = e_band + w_band from
tiltcrowd/wallshear (as in depth_correction.py); frames with exactly N = 2k-2 wall +1/2 in the second half;
rms(D_c) = circular sd of arg Z_N / N times the layer radius (cells); also the lock-corrected registry D = arg Z_N - tau (rad, twist channel).
K = 44, gamma1 = 60, chi = 0.5, lambda = 2. T_eff = (rms_meas / rms_pred(T = 1))^2. -> registry_stiffness.txt"""
import json, os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle
RT = "/home/baloewe/Physics/Overcharged/runs/cavity/server_main/main_act0.07"; K, G1, CHI, LAM = 44.0, 60.0, 0.5, 2.0
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/tiltcrowd"); import glob, defect_flow as DF, wallshear as WS
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"
SETS = {"passive (imposed vortex)": sorted(glob.glob(ROOT + "/runs/swirl_crowd/k0[468]_s*") + glob.glob(ROOT + "/runs/swirl_crowd/A03_k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A0.0[012]*")),
        "active": sorted(os.path.dirname(p) for p in glob.glob(RT + "/rho*/k0[3-9]/s0[1-3]/flowfield.dat.gz") + glob.glob(RT + "/rho*/k10/s0[1-3]/flowfield.dat.gz"))}
def meta(run):
    try:
        a = DF.one(run); w = WS.one(run); return run, a["k"], a["a"], a["depth"], w["e_band"] + w["w_band"]
    except Exception:
        return None
L = ["run | k | spacing a | wall shear s | l_E | k_D per defect | N | rms D_c meas (cells) | pred (T=1) | T_eff | sd tau (deg) | sd lock-corrected D (rad) | frames"]
res = []
for setname, runs in SETS.items():
  with Pool(6) as pool: M = [m for m in pool.map(meta, runs) if m]
  L.append(f"-- {setname}: {len(M)} runs")
  for rd, k, a, d, s in M:
    k = int(round(k)); bc = json.load(open(rd + "/input.json"))["BC"][0]; R = bc["R"]; c = np.array(bc["Q"][:2]); N = 2 * k - 2
    t, fr = read_with_angle(rd, c, 0.5); zs, taus, rs = [], [], []
    for f in fr:
        q = f[f[:, 2] > 0]; ph = np.arctan2(q[:, 1], q[:, 0]); r = np.hypot(q[:, 0], q[:, 1]); m = R - r < 4
        if m.sum() != N: continue
        zs.append(np.exp(1j * N * ph[m]).mean()); taus.append(np.angle(np.mean(np.exp(1j * (q[m, 3] - ph[m]))))); rs.append(r[m].mean())
    if len(zs) < 50 or not np.isfinite(s) or s <= 0: continue
    zs, taus = np.array(zs), np.array(taus); az = np.angle(zs); csd = lambda x: np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * x)))))
    lE = np.sqrt(K / (CHI * G1 * s * np.sqrt(LAM ** 2 - 1))); kD = K * np.pi ** 2 / (a * lE)
    meas = csd(az) / N * np.mean(rs); pred = np.sqrt(1.0 / (N * kD)); Te = (meas / pred) ** 2
    res.append((setname, k, a, Te)); L.append(f"{os.path.relpath(rd, ROOT + '/runs'):40s} | {k:2d} | {a:5.1f} | {s:.3f} | {lE:4.2f} | {kD:5.2f} | {N:2d} | {meas:.3f} | {pred:.3f} | {Te:6.1f} | "
                                     f"{np.degrees(csd(taus)):5.1f} | {csd(az - taus):.3f} | {len(zs)}")
for setname in SETS:
    R_ = [(k, a, Te) for sn, k, a, Te in res if sn == setname]
    if not R_: continue
    A_ = np.array(R_); L.append(f"\n{setname}: T_eff median {np.median(A_[:,2]):.1f} (IQR {np.percentile(A_[:,2],25):.1f}-{np.percentile(A_[:,2],75):.1f}); by k: " +
        ", ".join(f"{int(k)}: {np.median(A_[A_[:,0]==k,2]):.1f}" for k in sorted(set(A_[:,0]))) + f"; corr(log T_eff, log spacing) {np.corrcoef(np.log(A_[:,2]), np.log(A_[:,1]))[0,1]:+.2f}")
open("registry_stiffness.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L[-12:]))
