"""Director rotation vs depth in the one-way model: alpha(r) = director angle - local tangent (phi + 90 deg), as a nematic mean (mod 180), azimuthally
averaged, and the same at the defect angles vs between defects. Also S(r). Runs lam 2 and 4, U = 1, plus no flow; saves profile_dir.txt."""
import numpy as np, shear_tilt as M
np.seterr(all="ignore")
def state(lam, U, ring=12, T0=300, T1=1500, dt=0.04):
    zs = list((M.R - 2.6) * np.exp(2j * np.pi * (np.arange(ring) + 0.5) / ring)); th2, _, _ = M.fields(M.X, M.Y, zs, [0.5] * ring, M.R)
    q1 = np.where(M.inside, 0.5 * M.S0 * np.cos(th2), M.qb1); q2 = np.where(M.inside, 0.5 * M.S0 * np.sin(th2), M.qb2)
    for _ in range(int(T0 / dt)): q1, q2 = M.step(q1, q2, dt, 0, None)
    if U:
        fl = M.setup_flow(U, "noslip")
        for _ in range(int(T1 / dt)): q1, q2 = M.step(q1, q2, dt, lam * M.S0, fl)
    return q1, q2
L = []
for lam, U in ((2.0, 0.0), (2.0, 1.0), (4.0, 1.0), (1.2, 1.0)):
    q1, q2 = state(lam, U); th2 = np.arctan2(q2, q1); S = 2 * np.hypot(q1, q2); rel2 = th2 - 2 * (M.PH + np.pi / 2)   # 2 alpha
    o = M.observe(q1, q2, 1); L.append(f"lam {lam} U {U}: wall {o['nwall']} tilt {o.get('tilt', np.nan):.1f} coh {o.get('coh', np.nan):.2f}")
    for d in (0.25, 0.75, 1.5, 2.5, 3.5, 5, 7, 10, 14, 18, 24):
        m = M.inside & (abs((M.R - M.Rr) - d) < 0.3); z = np.mean(np.exp(1j * rel2[m]) * S[m])
        L.append(f"   depth {d:5.2f}: alpha = {np.degrees(np.angle(z)) / 2:+6.1f} deg (order |<S e^2ia>| {abs(z):.2f}), S {S[m].mean():.2f}")
print("\n".join(L)); open("profile_dir.txt", "w").write("\n".join(L) + "\n")
