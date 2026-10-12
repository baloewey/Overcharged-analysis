#!/usr/bin/env python3
"""E(tau) beyond first order in K1 - K3, with rigid cores (fixes the soft-core artefact of relax_frank.py, whose K1 = K3 check failed).
Each defect core is a disk of radius a in which the correction psi is a constant c_j (a rigid rotation of that defect's texture); the energy is
counted only outside the cores (hard mask). Variants: (A) c_j = 0: tilt imposed by the registry (tau); (B) c_j free: cores reorient, the
relaxed tilt is tau + 2 c_j. Same field, energy, grid and configuration as relax_frank.py. Checks: K1 = K3 -> flat E(tau), psi = 0.
Writes relax_frank2.txt and E_tau2_*.txt."""
import os, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
H = os.path.dirname(os.path.abspath(__file__))
ns = {"__file__": os.path.join(H, "relax_frank.py")}; exec(open(os.path.join(H, "relax_frank.py")).read().split("def energy_fn")[0], ns)
fields, config = ns["fields"], None
exec(open(os.path.join(H, "relax_frank.py")).read().split("def config")[1].split("def run")[0].join(["def config", ""]) if False else "", ns)
src = open(os.path.join(H, "relax_frank.py")).read(); exec("def config" + src.split("def config")[1].split("def run")[0], ns); config = ns["config"]
def run(k=6, R=31.830989, K1=42.0, K3=23.1, h=0.5, a=1.5, ntau=16, bulk="point", free_cores=False):
    x = np.arange(-R, R + h / 2, h); X, Y = np.meshgrid(x, x, indexing="ij"); n = len(x)
    I = sp.identity(n, format="csr"); D1 = sp.diags([-1, 1], [-1, 1], shape=(n, n)) / (2 * h)
    Dx = sp.kron(D1, I, format="csr"); Dy = sp.kron(I, D1, format="csr"); Kb, De = 0.5 * (K1 + K3), 0.5 * (K1 - K3)
    out = []
    for tau in np.linspace(0, 2 * np.pi, ntau, endpoint=False):
        zs, qs = config(k, R, tau, bulk); Z = X + 1j * Y
        inside = np.abs(Z) < R - 0.5 * h
        core_id = -np.ones(X.shape, int)
        for j, zi in enumerate(zs): core_id[(np.abs(Z - zi) < a) & inside] = j
        free = inside & (core_id < 0); nF = free.sum(); nc = len(zs)
        # psi_full = P psi_free + Q c   (Q: core indicator columns)
        P = sp.csr_matrix((np.ones(nF), (np.flatnonzero(free.ravel()), np.arange(nF))), shape=(X.size, nF))
        rows = np.flatnonzero((core_id >= 0).ravel()); Qm = sp.csr_matrix((np.ones(len(rows)), (rows, core_id.ravel()[rows])), shape=(X.size, nc))
        Th2, Gx, Gy = fields(X, Y, zs, qs, R); Th2, Gx, Gy = Th2.ravel(), Gx.ravel(), Gy.ravel()
        w = free.ravel().astype(float)                                   # energy only outside the cores
        def f(v):
            psi_f = v[:nF]; c = v[nF:] if free_cores else np.zeros(nc)
            pf = P @ psi_f + Qm @ c; gx = Gx + Dx @ pf; gy = Gy + Dy @ pf; C = np.cos(Th2 + 2 * pf); S = np.sin(Th2 + 2 * pf)
            dens = Kb * (gx * gx + gy * gy) + De * (-C * (gx * gx - gy * gy) - 2 * S * gx * gy); F = 0.5 * h * h * np.sum(w * dens)
            dgx = h * h * w * (Kb * gx + De * (-C * gx - S * gy)); dgy = h * h * w * (Kb * gy + De * (C * gy - S * gx))
            dpsi = h * h * w * De * (S * (gx * gx - gy * gy) - 2 * C * gx * gy)
            gfull = Dx.T @ dgx + Dy.T @ dgy + dpsi
            g = P.T @ gfull
            if free_cores: g = np.concatenate([g, Qm.T @ gfull])
            return F, g
        v0 = np.zeros(nF + (nc if free_cores else 0)); F0, _ = f(v0)
        res = minimize(f, v0, jac=True, method="L-BFGS-B", options=dict(maxiter=20000, gtol=1e-8))
        cw = res.x[nF:nF + 2 * k - 2] if free_cores else np.zeros(2 * k - 2)
        out.append((tau, F0, res.fun, tau + 2 * np.mean(cw), res.nit, res.success))
    return np.array(out, float)
def harm(E):
    c = np.fft.rfft(E) / len(E) * 2; return c[1].real, c[2].real, c[3].real
if __name__ == "__main__":
    L = []
    def say(s): print(s, flush=True); L.append(s)
    cases = [("check K1=K3, rigid cores fixed", dict(K1=32.5, K3=32.5)), ("check K1=K3, cores free", dict(K1=32.5, K3=32.5, free_cores=True)),
             ("A: cores fixed (tilt imposed), a 1.5", dict()), ("A: cores fixed, a 2.5", dict(a=2.5)), ("A: cores fixed, a 1.5, h 0.35", dict(h=0.35)),
             ("A: cores fixed, a 1.5, bulk pair", dict(bulk="pair")),
             ("B: cores free, a 1.5", dict(free_cores=True)), ("B: cores free, a 2.5", dict(a=2.5, free_cores=True))]
    for lab, kw in cases:
        t0 = time.time(); r = run(**kw); t, F0, F, trel = r.T[:4]
        a0, a1 = harm(F0), harm(F)
        say(f"{lab:40s}: unrelaxed cos tau {a0[0]:+7.2f} cos 2tau {a0[1]:+6.2f} | relaxed cos tau {a1[0]:+7.2f} cos 2tau {a1[1]:+6.2f} cos 3tau {a1[2]:+6.2f} kBT; "
            f"relaxation {np.mean(F0 - F):6.1f} kBT; relaxed tilt - tau rms {np.degrees(np.sqrt(np.mean(np.angle(np.exp(1j * (trel - t))) ** 2))):4.1f} deg; "
            f"conv {int(r[:, 5].sum())}/{len(r)}; {time.time() - t0:.0f} s")
        np.savetxt(os.path.join(H, "E_tau2_" + "".join(ch if ch.isalnum() else "_" for ch in lab) + ".txt"), r, header="tau F0 F tilt_rel nit ok")
    open(os.path.join(H, "relax_frank2.txt"), "w").write("\n".join(L) + "\n")
