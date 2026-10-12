#!/usr/bin/env python3
"""
Elastic (one-constant) defect interactions in an ellipse with strong anchoring whose winding is uniform in arclength (eMode 2).

Field: psi conjugate to the director angle; defects of charge q are sources, grad^2 psi = 2 pi sum q delta. Strong anchoring with
uniform winding per arclength fixes d(theta)/ds = 2 pi k / P on the wall, i.e. a Neumann condition with uniform flux for psi.
Neumann Green's function with uniform flux: G(z, w) = Log|z - w| + h_w(z), grad^2 h = 0, dh/dn = 2 pi / P - d/dn Log|z - w| on the wall
(compatible: the flux of Log|z - w| through the wall is 2 pi). Energy (units of pi K):
    E = - sum_{i != j} q_i q_j G(z_i, z_j) - sum_i q_i^2 h_{z_i}(z_i)        (core energies dropped)
For the disk this is exactly Log|z-w| + Log|R^2 - z conj(w)| + const (images of the same sign).
h_w is computed with the method of fundamental solutions: h = sum_m c_m Log|z - y_m| + c_0, sources y_m on the wall scaled by 1.15,
Neumann data collocated at 2M wall points, least squares.
Outputs (theory/): ellipse_shielding_results.txt and figures.
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = __file__.rsplit("/", 1)[0]
LA = 22.21


class Domain:
    def __init__(self, A, B, M=240, scale=1.15):
        self.A, self.B = A, B
        t = np.linspace(0, 2 * np.pi, 4 * M, endpoint=False)
        x, y = A * np.cos(t), B * np.sin(t)
        self.zc = x + 1j * y                                       # collocation points
        nx, ny = x / A ** 2, y / B ** 2; nn = np.hypot(nx, ny); self.n = (nx + 1j * ny) / nn
        ts = np.linspace(0, 2 * np.pi, M, endpoint=False)
        self.ys = scale * (A * np.cos(ts) + 1j * B * np.sin(ts))   # MFS sources outside
        tt = np.linspace(0, 2 * np.pi, 20001)
        self.P = np.sum(np.abs(np.diff(A * np.cos(tt) + 1j * B * np.sin(tt))))
        d = self.zc[:, None] - self.ys[None, :]
        Mat = np.real(d * np.conj(self.n[:, None])) / np.abs(d) ** 2   # d/dn Log|z - y|
        self.pinv = np.linalg.pinv(Mat, rcond=1e-12)

    def h_coeffs(self, w):
        d = self.zc - w
        rhs = 2 * np.pi / self.P - np.real(d * np.conj(self.n)) / np.abs(d) ** 2
        return self.pinv @ rhs

    def h(self, c, z):
        z = np.atleast_1d(z)
        return np.sum(c[None, :] * np.log(np.abs(z[:, None] - self.ys[None, :])), axis=1)

    def G(self, z, w, c=None):
        c = self.h_coeffs(w) if c is None else c
        return np.log(np.abs(np.atleast_1d(z) - w)) + self.h(c, z)

    def inside(self, z):
        return (z.real / self.A) ** 2 + (z.imag / self.B) ** 2 < 1


def energy(dom, zs, qs):
    E = 0.0
    cs = [dom.h_coeffs(w) for w in zs]
    for i, (zi, qi) in enumerate(zip(zs, qs)):
        E -= qi * qi * dom.h(cs[i], zi)[0]
        for j in range(i + 1, len(zs)):
            E -= 2 * qi * qs[j] * dom.G(zi, zs[j], cs[j])[0]
    return E


def relax(dom, zs0, qs):
    def f(v):
        zs = v[0::2] + 1j * v[1::2]
        if not np.all(dom.inside(zs)): return 1e6
        return energy(dom, zs, qs)
    v0 = np.empty(2 * len(zs0)); v0[0::2], v0[1::2] = zs0.real, zs0.imag
    r = minimize(f, v0, method="L-BFGS-B", options=dict(maxiter=4000))
    return r.x[0::2] + 1j * r.x[1::2], r.fun


def ring_start(A, B, N, depth, phase=0.0):
    t = np.linspace(0, 2 * np.pi, 20001); z = A * np.cos(t) + 1j * B * np.sin(t)
    s = np.r_[0, np.cumsum(np.abs(np.diff(z)))]; P = s[-1]
    targets = (np.arange(N) + phase) * P / N
    tt = np.interp(targets, s, t)
    zw = A * np.cos(tt) + 1j * B * np.sin(tt); nrm = np.cos(tt) / A + 1j * np.sin(tt) / B; nrm /= np.abs(nrm)
    return zw - depth * nrm


def tip_ratio(A, B, zs, band=4.0):
    t = np.linspace(0, 2 * np.pi, 6000, endpoint=False); W = A * np.cos(t) + 1j * B * np.sin(t)
    seg = np.abs(np.roll(W, -1) - W); tip = np.abs(W.real) > A / np.sqrt(2)
    j = np.argmin(np.abs(zs[:, None] - W[None, :]), axis=1); d = np.abs(zs - W[j])
    w = d < band
    return (np.sum(w & tip[j]) / seg[tip].sum()) / (np.sum(w & ~tip[j]) / seg[~tip].sum())


out = open(f"{HERE}/ellipse_shielding_results.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.write(s + "\n")

# ---- check against the disk
R = 31.83
dom = Domain(R, R)
w = 12.0 + 7.0j; z = np.array([-5.0 + 20.0j, 25.0 - 3.0j])
exact = np.log(np.abs(z - w)) + np.log(np.abs(R ** 2 - z * np.conj(w)))
num = dom.G(z, w)
say("disk check: G_num - G_exact at two points (should be equal up to one constant):", np.round(num - exact, 5))
for N in (12, 14):
    zs, E = relax(dom, ring_start(R, R, N, 3.0), np.full(N, 0.5))
    say(f"disk N={N}: relaxed ring depth {R - np.mean(np.abs(zs)):.3f} cells; analytic {R * (1 - ((N - 1) / (3 * N - 1)) ** (1 / (2 * N))):.3f}")

# ---- passive equilibrium in ellipses (area 2500, k = 6: N = 12 +1/2) and tip/flank ratio
area = 2500.0; N = 12; qs = np.full(N, 0.5)
res = {}
for eps in (1.0, 1.2, 1.5, 2.0):
    A = np.sqrt(area * eps / np.pi); B = A / eps
    d = Domain(A, B)
    best = None
    for ph in (0.0, 0.5):
        zs, E = relax(d, ring_start(A, B, N, 2.0, ph), qs)
        if best is None or E < best[1]: best = (zs, E)
    zs = best[0]; res[eps] = (d, zs)
    tw = np.linspace(0, 2 * np.pi, 6000); W = A * np.cos(tw) + 1j * B * np.sin(tw)
    depth = np.min(np.abs(zs[:, None] - W[None, :]), axis=1)
    say(f"eps {eps}: passive equilibrium, tip/flank density ratio {tip_ratio(A, B, zs):.3f} (sqrt(eps) = {np.sqrt(eps):.3f}, conductor = {eps}); "
        f"depth {depth.min():.2f}-{depth.max():.2f} cells")

# ---- shielding: potential felt by a bulk test defect (same sign +1/2 and opposite sign -1/2) from the equilibrium layer + its own image
fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
for col, eps in enumerate((1.0, 1.2, 1.5, 2.0)):
    d, zs = res[eps]
    cs = [d.h_coeffs(w) for w in zs]
    xs = np.linspace(-d.A, d.A, 121); ys = np.linspace(-d.B, d.B, 81)
    X, Y = np.meshgrid(xs, ys); Z = (X + 1j * Y).ravel()
    ins = (Z.real / d.A) ** 2 + (Z.imag / d.B) ** 2 < 0.8 ** 2          # inside the 0.8 level
    Zi = Z[ins]
    layer = np.zeros(len(Zi))
    for w, c in zip(zs, cs): layer += 0.5 * (np.log(np.abs(Zi - w)) + d.h(c, Zi))
    for row, qt in enumerate((0.5, -0.5)):
        selfim = np.array([d.h(d.h_coeffs(z0), z0)[0] for z0 in Zi]) * qt * qt
        V = -2 * qt * layer - selfim
        Vm = np.full(len(Z), np.nan); Vm[ins] = V - np.nanmean(V)
        ax = axs[row, col]
        im = ax.imshow(Vm.reshape(X.shape), origin="lower", extent=[-d.A, d.A, -d.B, d.B], cmap="RdBu_r", vmin=-0.6, vmax=0.6)
        ax.plot(zs.real, zs.imag, "o", ms=4, color="k"); t = np.linspace(0, 2 * np.pi, 300)
        ax.plot(d.A * np.cos(t), d.B * np.sin(t), color="k", lw=0.8); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        rng = np.nanmax(V) - np.nanmin(V); lay = -2 * qt * layer; rng_l = lay.max() - lay.min()
        ax.set_title(f"ε = {eps}, test {'+' if qt > 0 else '−'}½: total range {rng:.2f}\n(layer part {rng_l:.2f}) πK", fontsize=8)
        say(f"eps {eps}, test charge {qt:+}: energy range inside the 0.8 level: total {rng:.3f} pi K, from the layer alone {rng_l:.3f} pi K")
fig.colorbar(im, ax=axs, fraction=0.02, label="energy of the test defect (πK), mean removed")
fig.suptitle("Energy of a bulk test defect inside the 0.8 level, from 12 passive +½ at their equilibrium positions (black) with strong anchoring "
             "(uniform winding per arclength) and the test defect's own image; area 2500", fontsize=9)
fig.savefig(f"{HERE}/ellipse_interior_potential.png", dpi=150, bbox_inches="tight")
out.close()
