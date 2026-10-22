#!/usr/bin/env python3
"""Which wall defect of an ellipse is least bound, for one- and two-vortex tilt patterns? (user's question, 2026-10-08)
Field: theta = theta_s + psi, theta_s = sum_i q_i arg(z - z_i) (no images), psi single-valued: fixed on the wall band (nodes just outside the ellipse) to the
value that makes theta equal the eMode-2 anchoring k 2 pi s/P (s arclength from (A, 0), ccw; continuous branch along the wall), free inside, constant c_j on the
rigid core of defect j (radius a). Energy: anisotropic Frank (K1 = 42, K3 = 23.1; also one-constant 32.5), inside the ellipse, outside the cores; L-BFGS.
Configuration: k = 6, area 2500, 10 wall +1/2 evenly spaced in arclength at depth 1.6, bulk +1/2 pair at (+-0.45 A, 0).
Tilt patterns (tail tilt = tail direction - outward normal; tail = 2 theta_rest): natural (cores free); one vortex (+63 deg on every wall defect); two vortices
(+63 deg for x < 0, -63 deg for x > 0). Imposed tilts are set through c_j as the smallest rotation from the relaxed natural orientation of each core (avoids a pi branch jump between a core and its surroundings). Binding of wall defect j: move it inward along the normal by 1.5 cells
(same imposed tilt), re-relax: dE_j > 0; the smallest dE_j is the weakest-bound defect. Reported per eps and pattern: dE_j along the wall (by arclength
position, tips |x| > A/2 vs flanks), and the location of the weakest. Writes binding.txt and binding.png."""
import os, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); h = float(os.environ.get('GRID_H', 0.35)); a_core = 1.5; k = 6; DEPTH = 1.6; TILT = np.radians(63.0)
def wall_param(A, B, n=8000):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False); W = A * np.cos(t) + 1j * B * np.sin(t)
    seg = np.abs(np.roll(W, -1) - W); s = np.concatenate([[0], np.cumsum(seg)[:-1]]); P = seg.sum()
    T = (-A * np.sin(t) + 1j * B * np.cos(t)); T /= np.abs(T); Nn = -1j * T                    # outward normal
    return t, W, s, P, Nn
class Ellipse:
    def __init__(self, eps, off=(0.0, 0.0)):
        self.A = np.sqrt(2500 * eps / np.pi); self.B = self.A / eps; A, B = self.A, self.B
        x = np.arange(-A - 2, A + 2 + h / 2, h) + off[0]; y = np.arange(-B - 2, B + 2 + h / 2, h) + off[1]; X, Y = np.meshgrid(x, y, indexing="ij"); self.X, self.Y = X, Y; self.Z = X + 1j * Y
        nx, ny = len(x), len(y); Dx1 = sp.diags([-1, 1], [-1, 1], shape=(nx, nx)) / (2 * h); Dy1 = sp.diags([-1, 1], [-1, 1], shape=(ny, ny)) / (2 * h)
        self.Dx = sp.kron(Dx1, sp.identity(ny), format="csr"); self.Dy = sp.kron(sp.identity(nx), Dy1, format="csr")
        rho = np.sqrt((X / A) ** 2 + (Y / B) ** 2); self.inside = rho < 1
        self.t, self.W, self.s, self.P, self.Nn = wall_param(A, B)
        # wall band: outside nodes within 1.5 cells of the wall; nearest wall point index
        out = ~self.inside & (rho < 1 + 1.6 / min(A, B)); self.band = out
        idx = np.flatnonzero(out.ravel()); zb = self.Z.ravel()[idx]
        near = np.array([np.argmin(np.abs(self.W - zz)) for zz in zb]); self.band_idx, self.band_near = idx, near
    def site(self, sfrac, depth):
        i = int((sfrac % 1) * len(self.s)); i = np.searchsorted(self.s, (sfrac % 1) * self.P) % len(self.s)
        return self.W[i] - depth * self.Nn[i], self.Nn[i]
    def energy(self, zs, qs, tilt_target, K1=None, K3=None, free_cores=False, c_given=None):
        # defaults 42 / 23.1; env FRANK_K1 / FRANK_K3 override them (2026-10-22: one-constant K = 44 reruns of the escape energetics)
        if K1 is None: K1 = float(os.environ.get("FRANK_K1", 42.0))
        if K3 is None: K3 = float(os.environ.get("FRANK_K3", 23.1))
        """tilt_target: list (len = n_wall) of target tail tilts (rad) or None (natural); zs[:n_wall] are wall defects.
        c_given: optional dict j -> core constant (overrides tilt_target for those cores). Returns (energy, relaxed core constants)."""
        Z = self.Z; nw = len(tilt_target)
        e2 = np.ones(Z.shape, complex); gx = np.zeros(Z.shape); gy = np.zeros(Z.shape)
        for zi, q in zip(zs, qs):
            dz = Z - zi; dz = np.where(np.abs(dz) < 1e-9, 1e-9, dz); e2 *= (dz / np.abs(dz)) ** (2 * q); w = 1 / dz; gx += q * w.imag; gy += q * w.real
        Th2 = np.angle(e2).ravel(); Gx, Gy = gx.ravel(), gy.ravel()
        # boundary values of psi along the wall: anchoring theta = k 2 pi s / P minus theta_s, continuous branch
        th_s_wall = np.zeros(len(self.W))
        for zi, q in zip(zs, qs): th_s_wall += q * np.unwrap(np.angle(self.W - zi))
        g = k * 2 * np.pi * self.s / self.P - th_s_wall; g = np.unwrap(2 * g) / 2
        psi_b = np.zeros(Z.size); psi_b[self.band_idx] = g[self.band_near]
        core = -np.ones(Z.shape, int)
        for j, zi in enumerate(zs): core[(np.abs(Z - zi) < a_core) & self.inside] = j
        free = self.inside & (core < 0); nF = free.sum()
        Pm = sp.csr_matrix((np.ones(nF), (np.flatnonzero(free.ravel()), np.arange(nF))), shape=(Z.size, nF))
        rows = np.flatnonzero((core >= 0).ravel()); Qm = sp.csr_matrix((np.ones(len(rows)), (rows, core.ravel()[rows])), shape=(Z.size, len(zs)))
        # imposed core constants: tail tilt_j = 2 (theta_rest_s(z_j) + c_j) - normal angle  ->  c_j
        cfix = np.zeros(len(zs)); fixed = np.zeros(len(zs), bool)
        for j in range(nw):
            if c_given is not None and j in c_given: cfix[j] = c_given[j]; fixed[j] = True; continue
            if tilt_target[j] is None: continue
            zj = zs[j]; trs = sum(q * np.angle(zj - zi) for i, (zi, q) in enumerate(zip(zs, qs)) if i != j)
            nrm = np.angle(self.Nn[np.argmin(np.abs(self.W - zj))])
            cfix[j] = (tilt_target[j] + nrm) / 2 - trs; fixed[j] = True
        # cores not fixed: free variables (bulk pair natural, wall natural if None)
        freec = np.flatnonzero(~fixed); w = free.ravel().astype(float); Kb, De = 0.5 * (K1 + K3), 0.5 * (K1 - K3)
        base = psi_b + Qm @ cfix
        def f(v):
            pf = base + Pm @ v[:nF] + (Qm[:, freec] @ v[nF:] if len(freec) else 0)
            gxx = Gx + self.Dx @ pf; gyy = Gy + self.Dy @ pf; C = np.cos(Th2 + 2 * pf); S = np.sin(Th2 + 2 * pf)
            dens = Kb * (gxx ** 2 + gyy ** 2) + De * (-C * (gxx ** 2 - gyy ** 2) - 2 * S * gxx * gyy); F = 0.5 * h * h * np.sum(w * dens)
            dgx = h * h * w * (Kb * gxx + De * (-C * gxx - S * gyy)); dgy = h * h * w * (Kb * gyy + De * (C * gyy - S * gxx)); dps = h * h * w * De * (S * (gxx ** 2 - gyy ** 2) - 2 * C * gxx * gyy)
            gf = self.Dx.T @ dgx + self.Dy.T @ dgy + dps
            return F, np.concatenate([Pm.T @ gf, (Qm[:, freec].T @ gf) if len(freec) else []])
        v0 = np.zeros(nF + len(freec)); res = minimize(f, v0, jac=True, method="L-BFGS-B", options=dict(maxiter=30000, gtol=1e-8))
        call = cfix.copy(); call[freec] = res.x[nF:]
        self.last_trs = [sum(q * np.angle(zs[j] - zi) for i, (zi, q) in enumerate(zip(zs, qs)) if i != j) for j in range(nw)]
        self.last_nrm = [np.angle(self.Nn[np.argmin(np.abs(self.W - zs[j]))]) for j in range(nw)]
        return res.fun, call
def config(E, inward=None, delta=1.5):
    zs, normals, sfr = [], [], []
    for j in range(2 * k - 2):
        sf = (j + 0.5) / (2 * k - 2); z, n = E.site(sf, DEPTH)
        if inward == j: z = z - delta * n
        zs.append(z); normals.append(n); sfr.append(sf)
    zs += [0.45 * E.A + 0j, -0.45 * E.A + 0j]
    return zs, [0.5] * (2 * k), np.array(zs[:2 * k - 2])
if __name__ == "__main__":
    L = []
    def say(s_): print(s_, flush=True); L.append(s_)
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    for ie, eps in enumerate((1.0001, 1.5, 2.0)):
        offs = [(0.0, 0.0), (0.5 * h, 0.0), (0.0, 0.5 * h), (0.5 * h, 0.5 * h)]
        Es = [Ellipse(eps, o) for o in offs]; E = Es[0]; t0 = time.time(); _, _, zw = config(E)
        pats = {"natural": [None] * 10, "one vortex": [TILT] * 10, "two vortices": [TILT if z.real < 0 else -TILT for z in zw]}
        def imposed(Eo, zs, qs, pat):
            En, cn = Eo.energy(zs, qs, [None] * 10)
            if all(p_ is None for p_ in pat): return En
            tilt_nat = [2 * (Eo.last_trs[j] + cn[j]) - Eo.last_nrm[j] for j in range(10)]
            cg = {j: cn[j] + np.angle(np.exp(1j * (pat[j] - tilt_nat[j]))) / 2 for j in range(10)}
            return Eo.energy(zs, qs, [None] * 10, c_given=cg)[0]
        for pname, pat in pats.items():
            per_off = []
            for Eo in Es:
                zs, qs, _ = config(Eo); E0 = imposed(Eo, zs, qs, pat); dE = []
                for j in range(10):
                    zs_j, qs_j, _ = config(Eo, inward=j); dE.append(imposed(Eo, zs_j, qs_j, pat) - E0)
                per_off.append(dE)
            per_off = np.array(per_off); dE = per_off.mean(0); err = per_off.std(0) / np.sqrt(len(offs))
            tips = np.abs(zw.real) > E.A / 2; jmin = int(np.argmin(dE))
            say(f"eps {eps:.1f} {pname:12s}: binding dE_j (kBT; mean over 4 grid offsets +- s.e.) by site (x/A): " + " ".join(f"{z.real / E.A:+.2f}:{d:+.1f}+-{e:.1f}" for z, d, e in zip(zw, dE, err)) +
                f" | tips {dE[tips].mean() if tips.any() else np.nan:+.1f}, flanks {dE[~tips].mean():+.1f}; site spread (sd) {dE.std():.1f}; weakest at x/A = {zw[jmin].real / E.A:+.2f} ({'tip' if tips[jmin] else 'flank'})")
            ax[ie].errorbar(np.angle(zw), dE, yerr=err, fmt="o-", ms=3, label=pname)
        ax[ie].set_title(f"$\\varepsilon$ = {eps:.1f}", fontsize=10); ax[ie].set_xlabel("polar angle of the site (rad; 0, $\\pm\\pi$ = tips)"); ax[ie].set_ylabel("binding $\\Delta E_j$ ($k_BT$)")
        say(f"   ({time.time() - t0:.0f} s)")
    ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(H, "binding.png"), dpi=150); open(os.path.join(H, "binding.txt"), "w").write("\n".join(L) + "\n")
