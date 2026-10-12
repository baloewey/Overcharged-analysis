#!/usr/bin/env python3
"""Is a uniformly tilted wall layer cheaper with two defects ejected? (user's mechanism, 2026-10-08: vortex -> flow -> flow tilts the wall defects
-> the 2k ring becomes too expensive -> two +1/2 leave into the bulk.)
Structures (same 2k defects of charge +1/2 in both, so core energies cancel):
  A: 2k +1/2 on a ring at the passive depth-law radius;
  B: 2k - 2 +1/2 on a ring at the active depth-law radius (bulk charge 1), plus two +1/2 at +-s on a diameter.
Imposed uniform tilt tau of every wall defect (tail angle - outward normal): each core is a rigid disk (radius a) whose director is the one-constant
field plus a constant c_j = (tau - tau_j^0)/2, tau_j^0 = the natural tail tilt of defect j in the one-constant field (the bulk pair keeps its natural
orientation, c = 0); B's ring is first shifted by the registry alpha = tau/(2k - 2) (the lock), so its c_j are ~0. Outside the cores the correction
psi (psi = 0 on the wall, = c_j on core j) is relaxed under the Frank energy, one-constant (K1 = K3 = 32.5) or anisotropic (K1 = 42, K3 = 23.1).
Validation: tau = 0, one-constant, unrelaxed vs the analytic image energy difference.
Output: E_A(tau), E_B(tau), Delta E = E_B - E_A (kBT), for tau = 0..90 deg, and s = 0.3, 0.45, 0.6 R. Writes tilted_structures.txt and .png."""
import os, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
def fields(X, Y, zs, qs, R):
    z = X + 1j * Y; e2 = np.ones(z.shape, complex); gx = np.zeros(z.shape); gy = np.zeros(z.shape)
    for zi, q in zip(zs, qs):
        a = z - zi; a = np.where(np.abs(a) < 1e-9, 1e-9, a); b = 1 - z * np.conj(zi) / R ** 2
        e2 *= ((a / np.abs(a)) * (b / np.abs(b))) ** (2 * q); wv = 1 / a - (np.conj(zi) / R ** 2) / b
        gx += q * wv.imag; gy += q * wv.real
    return np.angle(e2), gx, gy
def natural_tilts(zs, qs, R, n_wall, eps=1e-3):
    out = []
    for j in range(n_wall):
        zj = zs[j]; p = zj + eps * zj / abs(zj)
        e2, _, _ = fields(np.array([p.real]), np.array([p.imag]), zs, qs, R)
        head2 = e2[0] - np.angle(p - zj)                      # 2 theta_rest (tail direction)
        out.append(np.angle(np.exp(1j * (head2 - np.angle(zj)))))
    return np.array(out)
def analytic_E(zs, qs, R):                                        # one-constant image energy, units of pi K
    E = 0.0; n = len(zs)
    for i in range(n):
        for j in range(n):
            if i != j: E -= qs[i] * qs[j] * np.log(abs(zs[i] - zs[j]))
            E -= qs[i] * qs[j] * np.log(abs(R ** 2 - zs[i] * np.conj(zs[j])))
    return E
def energy(zs, qs, cfix, R, K1, K3, h=0.5, a=1.5, relax=True):
    x = np.arange(-R, R + h / 2, h); X, Y = np.meshgrid(x, x, indexing="ij"); n = len(x); Z = X + 1j * Y
    I = sp.identity(n, format="csr"); D1 = sp.diags([-1, 1], [-1, 1], shape=(n, n)) / (2 * h); Dx = sp.kron(D1, I, format="csr"); Dy = sp.kron(I, D1, format="csr")
    inside = np.abs(Z) < R - 0.5 * h; core = -np.ones(X.shape, int)
    for j, zi in enumerate(zs): core[(np.abs(Z - zi) < a) & inside] = j
    free = inside & (core < 0); nF = free.sum()
    P = sp.csr_matrix((np.ones(nF), (np.flatnonzero(free.ravel()), np.arange(nF))), shape=(X.size, nF))
    cfield = np.zeros(X.size); cr = np.flatnonzero((core >= 0).ravel()); cfield[cr] = np.array(cfix)[core.ravel()[cr]]
    Th2, Gx, Gy = fields(X, Y, zs, qs, R); Th2, Gx, Gy = Th2.ravel(), Gx.ravel(), Gy.ravel(); w = free.ravel().astype(float)
    Kb, De = 0.5 * (K1 + K3), 0.5 * (K1 - K3)
    def f(v):
        pf = P @ v + cfield; gx = Gx + Dx @ pf; gy = Gy + Dy @ pf; C = np.cos(Th2 + 2 * pf); S = np.sin(Th2 + 2 * pf)
        dens = Kb * (gx * gx + gy * gy) + De * (-C * (gx * gx - gy * gy) - 2 * S * gx * gy); F = 0.5 * h * h * np.sum(w * dens)
        dgx = h * h * w * (Kb * gx + De * (-C * gx - S * gy)); dgy = h * h * w * (Kb * gy + De * (C * gy - S * gx)); dpsi = h * h * w * De * (S * (gx * gx - gy * gy) - 2 * C * gx * gy)
        return F, P.T @ (Dx.T @ dgx + Dy.T @ dgy + dpsi)
    v0 = np.zeros(nF); F0 = f(v0)[0]
    if not relax: return F0
    res = minimize(f, v0, jac=True, method="L-BFGS-B", options=dict(maxiter=20000, gtol=1e-8)); return res.fun
def structA(k, R):
    N = 2 * k; rs = R * ((N - 1) / (3 * N - 1)) ** (1 / (2 * N)) - 0.61
    return list(rs * np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)), [0.5] * N, N
def structB(k, R, s, tau):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; al = tau / N
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [s + 0j, -s + 0j], [0.5] * (N + 2), N
def cores_for(zs, qs, R, nw, tau):
    t0 = natural_tilts(zs, qs, R, nw); c = list(np.angle(np.exp(1j * (tau - t0))) / 2) + [0.0] * (len(zs) - nw); return c, t0
if __name__ == "__main__":
    L = []
    def say(s_): print(s_, flush=True); L.append(s_)
    k, R = 6, 31.830989
    # validation
    zA, qA, nA = structA(k, R); zB, qB, nB = structB(k, R, 0.45 * R, 0.0)
    eA = energy(zA, qA, [0.0] * len(zA), R, 32.5, 32.5, relax=False); eB = energy(zB, qB, [0.0] * len(zB), R, 32.5, 32.5, relax=False)
    say(f"validation (tau 0 natural cores, one-constant, unrelaxed): numerical E_B - E_A = {eB - eA:+.1f} kBT; analytic pi K dE = {np.pi * 32.5 * (analytic_E(zB, qB, R) - analytic_E(zA, qA, R)):+.1f} kBT "
        f"(the analytic value includes the core self-energies' position dependence; differences reflect the core cutoff)")
    say(f"natural tail tilts: ring A (2k = {nA}) {np.round(np.degrees(natural_tilts(zA, qA, R, nA)), 0)}; ring B (2k-2 = {nB}) {np.round(np.degrees(natural_tilts(zB, qB, R, nB)), 0)}")
    taus = np.radians([0, 15, 30, 45, 60, 75, 90]); res = {}
    for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
        for sfac in (0.3, 0.45, 0.6):
            EA, EB = [], []
            t0 = time.time()
            for tau in taus:
                zA, qA, nA = structA(k, R); cA, _ = cores_for(zA, qA, R, nA, tau)
                zB, qB, nB = structB(k, R, sfac * R, tau); cB, _ = cores_for(zB, qB, R, nB, tau)
                EA.append(energy(zA, qA, cA, R, K1, K3)); EB.append(energy(zB, qB, cB, R, K1, K3))
            EA, EB = np.array(EA), np.array(EB); res[(lab, sfac)] = (EA, EB)
            say(f"{lab:12s} s = {sfac:.2f} R: E_A(tau) - E_A(0) = " + " ".join(f"{np.degrees(t):.0f}:{e:+.0f}" for t, e in zip(taus, EA - EA[0])) +
                " | E_B(tau) - E_B(0) = " + " ".join(f"{np.degrees(t):.0f}:{e:+.0f}" for t, e in zip(taus, EB - EB[0])) +
                " | Delta E = E_B - E_A = " + " ".join(f"{np.degrees(t):.0f}:{e:+.0f}" for t, e in zip(taus, EB - EA)) + f" kBT  ({time.time() - t0:.0f} s)")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for i, lab in enumerate(("one-constant", "anisotropic")):
        for sfac in (0.3, 0.45, 0.6):
            EA, EB = res[(lab, sfac)]; ax[i].plot(np.degrees(taus), EB - EA, "o-", label=f"bulk pair at $\\pm{sfac}R$")
        ax[i].axhline(0, color="k", lw=0.6); ax[i].set_xlabel("imposed uniform tilt (deg)"); ax[i].set_ylabel(r"$E_B-E_A$ ($k_BT$)"); ax[i].set_title(f"({'ab'[i]}) {lab}", fontsize=10); ax[i].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(H, "tilted_structures.png"), dpi=150); open(os.path.join(H, "tilted_structures.txt"), "w").write("\n".join(L) + "\n")
