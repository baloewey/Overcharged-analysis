#!/usr/bin/env python3
"""Cassini-oval geometry for the elastic point-defect model (binding.Ellipse.energy is reused unchanged). User, 2026-10-12.
Wall: (x^2+y^2)^2 - 2a^2(x^2-y^2) = b^4 - a^4, i.e. d1 d2 = b^2 with foci (+-a, 0) (solver "shape": 1, "cassiniA" = a, "R" = b). Anchoring as eMode 3: director angle
k 2 pi s/P along the arclength s from (x_max, 0), counter-clockwise (same rule as eMode 2 for ellipses). Wall points from the polar parametrisation (as wall_fourier.wall_polyline).
CassiniAB(a, b) provides the attributes that binding.Ellipse.energy uses: Z, Dx, Dy, inside, W (complex), s, P, Nn (outward unit normals, complex), band_idx, band_near;
plus contains(z) for point tests and site(sfrac, depth) (inherited)."""
import os, sys, numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
TH = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path.insert(0, TH + "/ellipse_binding"); os.environ.setdefault("GRID_H", "0.5")
import binding as bd  # noqa
def cassini_wall(a, b, n=8000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False); r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    W = np.sqrt(r2) * np.exp(1j * th); seg = np.abs(np.roll(W, -1) - W); s = np.concatenate([[0], np.cumsum(seg)[:-1]]); P = seg.sum()
    # outward normal = grad of F = |z-a|^2 |z+a|^2 (direction only)
    x, y = W.real, W.imag; d1 = (x - a) ** 2 + y ** 2; d2 = (x + a) ** 2 + y ** 2
    gx = 2 * (x - a) * d2 + 2 * (x + a) * d1; gy = 2 * y * d2 + 2 * y * d1; Nn = (gx + 1j * gy) / np.hypot(gx, gy)
    return W, s, P, Nn
class CassiniAB(bd.Ellipse):
    def __init__(self, a, b):
        self.a, self.b = a, b; h = bd.h
        W, s, P, Nn = cassini_wall(a, b); self.W, self.s, self.P, self.Nn = W, s, P, Nn; self.t = np.angle(W)
        xm, ym = np.abs(W.real).max(), np.abs(W.imag).max(); self.A, self.B = xm, ym
        x = np.arange(-xm - 2, xm + 2 + h / 2, h); y = np.arange(-ym - 2, ym + 2 + h / 2, h); X, Y = np.meshgrid(x, y, indexing="ij"); self.X, self.Y = X, Y; self.Z = X + 1j * Y
        nx, ny = len(x), len(y); Dx1 = sp.diags([-1, 1], [-1, 1], shape=(nx, nx)) / (2 * h); Dy1 = sp.diags([-1, 1], [-1, 1], shape=(ny, ny)) / (2 * h)
        self.Dx = sp.kron(Dx1, sp.identity(ny), format="csr"); self.Dy = sp.kron(sp.identity(nx), Dy1, format="csr")
        self.inside = self.contains(self.Z)
        tree = cKDTree(np.c_[W.real, W.imag]); dist, near = tree.query(np.c_[X.ravel(), Y.ravel()]); dist = dist.reshape(X.shape); near = near.reshape(X.shape)
        out = ~self.inside & (dist < 1.6); self.band = out; idx = np.flatnonzero(out.ravel()); self.band_idx, self.band_near = idx, near.ravel()[idx]
    def contains(self, z):
        z = np.asarray(z); return np.abs(z - self.a) * np.abs(z + self.a) < self.b ** 2
if __name__ == "__main__":
    # validation: a = 0 is a circle of radius b; compare with binding.Ellipse(1.0) (area 2500) for the same defect configuration
    bd.k = 6; E0 = bd.Ellipse(1.0); R = E0.A; C = CassiniAB(0.0, R)
    zs = [C.site((j + 0.5) / 10, 1.6)[0] for j in range(10)] + [0.45 * R, -0.45 * R]; qs = [0.5] * 12
    e_c = C.energy(zs, qs, [None] * 12)[0]; e_e = E0.energy([E0.site((j + 0.5) / 10, 1.6)[0] for j in range(10)] + [0.45 * R, -0.45 * R], qs, [None] * 12)[0]
    print(f"circle via Cassini a=0: {e_c:.1f} kBT; via binding.Ellipse: {e_e:.1f} kBT; wall perimeter {C.P:.2f} vs {E0.P:.2f}")
    # an elongated oval: inside fraction and perimeter
    C2 = CassiniAB(22.6186, 24.8805); print(f"oval a 22.62 b 24.88: x_max {C2.A:.2f}, y_max {C2.B:.2f}, area {C2.inside.sum() * bd.h ** 2:.0f}, perimeter {C2.P:.1f}")
