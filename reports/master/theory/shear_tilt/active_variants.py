#!/usr/bin/env python3
"""Continuum tests of the MPCD crowding drop (user, 2026-10-19): is the drop due to noise or active turbulence that the deterministic, laminar continuum lacks?
Active Q-tensor + backflow (one constant, active_bf.run), R 31.83, lambda 2, k 4 and 16 (2k-2 ring + bulk pair seeded).
  noise: zeta 0.3 with additive noise sqrt(2 D dt) xi on q1, q2 inside the cavity after every step, D = 0.002 or 0.01;
  activity: zeta 3 (dt 0.005) and 10 (dt 0.002), no noise.
Records: wall count, coherence, layer tilt (last observations), circulation, bulk +-1/2 counts, and the bulk spiral order (azimuthal <S e^{2i(theta-phi)}>, 0.3-0.6 R).
Usage: active_variants.py MODE k PARAM  (MODE noise|act) -> var_<MODE>_<PARAM>_k<k>.json"""
import sys, json, numpy as np
import shear_tilt as M
mode, k, par = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
M.k = k; M.qb1 = 0.5 * M.S0 * np.cos(2 * k * M.PH); M.qb2 = 0.5 * M.S0 * np.sin(2 * k * M.PH)
import active_bf as B
rng = np.random.default_rng(100 * k + int(par * 1000))
zeta, dt = (0.3, 0.02) if mode == "noise" else (par, 0.005 if par <= 3 else 0.002)
if mode == "noise":
    step0 = M.step
    def noisy(q1, q2, dt_, lam_c, fl):
        q1, q2 = step0(q1, q2, dt_, lam_c, fl)
        if fl is not None:
            a = np.sqrt(2 * par * dt_); q1 = q1 + a * rng.standard_normal(q1.shape) * M.inside; q2 = q2 + a * rng.standard_normal(q2.shape) * M.inside
        return q1, q2
    M.step = noisy
obs, (q1, q2, vx, vy) = B.run(2.0, zeta, 0.0, ring=2 * k - 2, T0=300, T1=3000, dt=dt, every=100, nflow=1 if mode == "noise" else 2)
last = obs[-10:]; zs = [np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]
r = np.hypot(M.X, M.Y); ph = np.arctan2(M.Y, M.X); m = M.inside & (r > 0.3 * M.R) & (r < 0.6 * M.R)
Z = np.mean(2 * np.hypot(q1, q2)[m] * np.exp(1j * (np.arctan2(q2, q1)[m] - 2 * ph[m])))
out = dict(mode=mode, k=k, par=par, zeta=zeta, nwall=float(np.mean([o["nwall"] for o in last])), coh=float(np.nanmean([o.get("coh", np.nan) for o in last])),
           tilt=float(np.degrees(np.angle(np.mean(zs)))) if zs else None, circ=float(np.mean([o["circ"] for o in last])), nbulk_p=float(np.mean([o.get("nbulk_p", 0) for o in last])),
           nbulk_m=float(np.mean([o.get("nbulk_m", 0) for o in last])), spiral_order=float(abs(Z)), tilt_series=[o.get("tilt") for o in obs])
json.dump(out, open(f"var_{mode}_{par:g}_k{k:02d}.json", "w")); print(json.dumps({k_: v for k_, v in out.items() if k_ != "tilt_series"}))
