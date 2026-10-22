#!/usr/bin/env python3
"""Flow-induced rotation of the director near the wall in the continuum (user, 2026-10-22; companion of bandrot/band_rotation.py for MPCD).
Flat no-slip wall in simple shear (flatwall.py model, one constant, L = Gamma = 1, S0 0.876, lambda 2, one-way flow), one +1/2 per site spacing a_s.
After the sheared state has converged (3000 t.u.), dtheta = theta_LdG - theta_static, where theta_static is the exact static one-constant field
of a +1/2 row at the measured defect position (flatpin/flat_pinning.setup: row with images + harmonic correction, anchoring exact on the wall).
Profile against distance from the wall (0.5 bins), cells farther than 1.5 from the core; signed by the tilt (positive = rotated further in the
flow sense). The band estimate chi gamma1 s d^2/(4K) has chi = 1, gamma1/K = 1/(S0^2) ... reported in model units as gamma1 s d^2/(4K) with
gamma1/K = 1 (flatwall_depth convention). Shears s 0.03 and 0.06, a_s 8, 12, 20. -> flatwall_band.txt, flatwall_band.npz"""
import sys, numpy as np
from multiprocessing import Pool
import flatwall as F
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/flatpin"); import flat_pinning as FP
S0, h, c = F.S0, F.h, F.c; UB = np.arange(0, 10.01, 0.5)
def run(job):
    a_s, s, lam = job; H = 40.0; dt = 0.04; nx, ny = int(round(a_s / h)), int(round(H / h)) + 2
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 1) * h; X, Y = np.meshgrid(x, y, indexing="ij")
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2): q1[:, 0] = qw1; q2[:, 0] = qw2; q1[:, -1] = q1[:, -2]; q2[:, -1] = q2[:, -2]; return q1, q2
    lap = lambda f: (np.roll(f, 1, 0) + np.roll(f, -1, 0) - 2 * f) / h ** 2 + (np.roll(f, 1, 1) + np.roll(f, -1, 1) - 2 * f) / h ** 2
    gx = lambda f: (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
    vx = s * np.maximum(Y, 0); Exy = Oxy = 0.5 * s; lam_c = lam * S0
    def step(q1, q2, flow):
        tr = 2 * (q1 ** 2 + q2 ** 2); H1 = lap(q1) + q1 - c * tr * q1; H2 = lap(q2) + q2 - c * tr * q2
        if flow: d1 = -vx * gx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * gx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
        else: d1, d2 = H1, H2
        return bc(q1 + dt * d1, q2 + dt * d2)
    q1, q2 = bc(q1, q2)
    for _ in range(int(300 / dt)): q1, q2 = step(q1, q2, False)
    for _ in range(int(3000 / dt)): q1, q2 = step(q1, q2, True)
    m = F.measure(q1, q2, X, Y, a_s); 
    if "tilt" not in m: return job, None
    S = 2 * np.hypot(q1, q2); inn = Y > 0.3; import scipy.ndimage as nd
    lab, n = nd.label((S < 0.5 * S0) & inn, structure=np.ones((3, 3))); i = np.argmax(nd.sum(S0 - S, lab, range(1, n + 1))) + 1; mk = lab == i
    ang = np.angle(np.sum(np.exp(2j * np.pi * X[mk] / a_s) * (S0 - S[mk]))); cx = (ang % (2 * np.pi)) * a_s / (2 * np.pi); cy = float(np.sum(Y[mk] * (S0 - S[mk])) / np.sum(S0 - S[mk]))
    Xs, Ys, ths = FP.setup(a_s, cy, cx, hh=h, ymax=H)[:3]                     # static exact field on x = (i + 1/2) h, y = (j + 1/2) h
    tl = 0.5 * np.angle(q1 + 1j * q2); dth = np.full(X.shape, np.nan)
    for j in range(1, ny):
        yy = Y[0, j]
        if yy < h / 2 or yy > Ys[0, -1]: continue
        jj = np.searchsorted(Ys[0], yy); jj = min(max(jj, 1), Ys.shape[1] - 1); f = (yy - Ys[0, jj - 1]) / h
        e2 = (1 - f) * np.exp(2j * ths[:, jj - 1]) + f * np.exp(2j * ths[:, jj]); dth[:, j] = np.angle(np.exp(2j * tl[:, j]) / e2) / 2
    sg = np.sign(m["tilt"]) or 1.0; dxp = (X - cx + a_s / 2) % a_s - a_s / 2; far = np.hypot(dxp, Y - cy) > 1.5; sel = far & np.isfinite(dth) & (Y > 0)
    hst, _ = np.histogram(Y[sel], UB, weights=sg * dth[sel]); cnt, _ = np.histogram(Y[sel], UB)
    return job, dict(tilt=m["tilt"], depth=cy, prof=np.degrees(hst / np.maximum(cnt, 1)), est=np.degrees(s * cy ** 2 / 4))
if __name__ == "__main__":
    jobs = [(a_s, s, 2.0) for s in (0.03, 0.06) for a_s in (8.0, 12.0, 20.0)]
    with Pool(6) as p: res = p.map(run, jobs)
    uc = 0.5 * (UB[1:] + UB[:-1]); L = ["CONTINUUM FLAT WALL: dtheta = LdG - exact static field of the same defect (deg, + = further in the flow sense); band estimate s d^2/4 (model units, gamma1/K = 1)",
         "a_s | s | tilt | depth | band estimate | dtheta at u = 0.25, 0.75, 1.25, 1.75, 2.25 | 3.25, 4.75, 7.25"]; ix = lambda v: int(np.argmin(np.abs(uc - v)))
    for (a_s, s, lam), r in res:
        if r is None: L.append(f"{a_s} {s}: no defect"); continue
        p = r["prof"]; L.append(f"{a_s:4.0f} | {s:.2f} | {r['tilt']:5.1f} | {r['depth']:.2f} | {r['est']:+5.2f} | " + " ".join(f"{p[ix(v)]:+5.1f}" for v in (0.25, 0.75, 1.25, 1.75, 2.25)) + " | " + " ".join(f"{p[ix(v)]:+5.1f}" for v in (3.25, 4.75, 7.25)))
    np.savez("flatwall_band.npz", uc=uc, **{f"a{a_s:g}_s{s:g}": r["prof"] for (a_s, s, lam), r in res if r})
    open("flatwall_band.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
