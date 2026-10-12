#!/usr/bin/env python3
"""Flow imposed on a wall +1/2 by an approaching bulk +1/2. Disk R = 31.83, k = 6; director = one-constant texture with images (tilted_structures.fields) of the ring of 10
(depth 2.06) + a fixed bulk +1/2 at (-0.45 R, 0) + the approacher; active force f = -(zeta/eta) S/2 div(Q-hat) with zeta/eta = 0.5 (director-field Stokes fit 0.42-0.67) and S = 0.876;
no-slip Stokes (active_bf.stokes). The approacher's contribution = flow(approacher at P) - flow(approacher far at (0, 0.45 R)).
Reported at wall defect j (angle 0, depth 2.06): inward radial velocity, tangential velocity, and the change of the local shear rate dv_t/dr at j (finite difference over 1 cell),
vs approacher position (in front at r = 3..14, side). Compare: measured fluid speed at the layer 0.13-0.16, wall-defect inward drift 0.04-0.09 cells per 10 t.u. -> flow_kick.txt"""
import os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); TH = os.path.dirname(H)
sys.path.insert(0, os.path.join(TH, "ejection")); sys.path.insert(0, os.path.join(TH, "shear_tilt"))
from tilted_structures import fields  # noqa
import active_bf as B  # noqa  (disk grid, stokes)
M = B.M; R = M.R; N = 10; d0 = 2.06; rs = R - d0; zeta_eta = 0.5; S = 0.876
ring = [rs * np.exp(2j * np.pi * m / N) for m in range(N)]; fixed = -0.45 * R + 0j
def flow_for(app):
    zs = ring + [fixed, app]; th2, _, _ = fields(M.X, M.Y, zs, [0.5] * len(zs), R)
    q1 = np.where(M.inside, 0.5 * S * np.cos(th2), M.qb1); q2 = np.where(M.inside, 0.5 * S * np.sin(th2), M.qb2)
    fx = -zeta_eta * (B.dx(q1) + B.dy(q2)) * M.inside; fy = -zeta_eta * (B.dx(q2) - B.dy(q1)) * M.inside
    return B.stokes(fx, fy)
def at(vx, vy, z):
    from scipy.ndimage import map_coordinates
    ix = (z.real - M.x[0]) / M.h; iy = (z.imag - M.x[0]) / M.h
    return map_coordinates(vx, [[ix], [iy]], order=1)[0], map_coordinates(vy, [[ix], [iy]], order=1)[0]
def probe(vx, vy):
    zj = ring[0]; ux, uy = at(vx, vy, zj); er = np.array([1.0, 0.0]); et = np.array([0.0, 1.0])     # j at angle 0: radial = x, tangential = y
    inward = -(ux * er[0] + uy * er[1]); tang = ux * et[0] + uy * et[1]
    zo, zi = (rs + 0.5) + 0j, (rs - 0.5) + 0j; _, vto = at(vx, vy, zo); _, vti = at(vx, vy, zi); shear = (vto - vti) / 1.0
    return inward, tang, shear
if __name__ == "__main__":
    ref = probe(*flow_for(0.45j * R)); L = []
    L.append(f"reference (approacher far): at wall defect j: inward {ref[0]:+.4f}, tangential {ref[1]:+.4f}, shear dv_t/dr {ref[2]:+.4f}  (zeta/eta = {zeta_eta})")
    cases = [(f"in front, r = {r}", rs - r + 0j) for r in (3, 4, 5, 6, 8, 10, 14)] + [("side, depth 6, half a site along", (R - 6) * np.exp(1j * np.pi / N))]
    for name, z in cases:
        p = probe(*flow_for(z))
        L.append(f"approacher {name:34s}: extra inward {p[0] - ref[0]:+.4f}, extra tangential {p[1] - ref[1]:+.4f}, extra shear {p[2] - ref[2]:+.4f} (cells / t.u.)")
    print("\n".join(L)); open(os.path.join(H, "flow_kick.txt"), "w").write("\n".join(L) + "\n")
