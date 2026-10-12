#!/usr/bin/env python3
"""Passive wall lattice: exact Boltzmann sampling (Metropolis) of the full electrostatic energy in a disk, vs the harmonic (Hessian) result
and the measurement. E/(pi K) = - sum_{i!=j} q^2 ln|z_i - z_j| - sum_{i,j} q^2 ln|Rw^2 - z_i conj(z_j)|, q = 1/2, N = 2k (circle R = 31.83).
Variants: K in {32.5 (mean K1, K3), 23 (bend), 17.8 (K_eff fitted on the tangential rms)}; effective wall Rw = R or R - 0.61 (depth-law offset);
detector: exact positions, or rounded to the integer cell grid (defects.dat stores cell positions; grid offset from the run's centre Q).
Measured as in measure_passive.py: depth = R - r, internal tangential rms (arclength deviation about the site, collective removed).
Writes mc_exact.txt."""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
M = json.load(open(os.path.join(H, "fluct_measured.json")))
rng = np.random.default_rng(5)
R = 31.830989; Q = 33.5                         # Rscan_R31.83 centre (grid offset)
def energy_i(z, i, Rw):
    zi = z[i]; o = np.delete(z, i)
    return -0.25 * (2 * np.log(np.abs(zi - o)).sum() + 2 * np.log(np.abs(Rw ** 2 - zi * np.conj(o))).sum() + np.log(Rw ** 2 - abs(zi) ** 2))
def sample(N, K, Rw, sweeps=40000, burn=4000, every=10, step=0.6):
    beta = np.pi * K                            # energy unit pi K, kBT = 1
    depth0 = R * (1 - ((N - 1) / (3 * N - 1)) ** (1 / (2 * N)))
    z = (Rw - depth0) * np.exp(2j * np.pi * (np.arange(N) + 0.5) / N); out = []; acc = 0
    for sw in range(sweeps):
        for i in rng.permutation(N):
            old = z[i]; e0 = energy_i(z, i, Rw)
            z[i] = old + step * (rng.normal() + 1j * rng.normal())
            if abs(z[i]) >= Rw - 0.05: z[i] = old; continue
            e1 = energy_i(z, i, Rw)
            if rng.random() < np.exp(-beta * (e1 - e0)): acc += 1
            else: z[i] = old
        if sw >= burn and sw % every == 0: out.append(z.copy())
    return np.array(out), acc / (sweeps * N)
def measure(Z, N, rounded):
    if rounded:
        Z = (np.round(Z.real + Q) - Q) + 1j * (np.round(Z.imag + Q) - Q)
    r = np.abs(Z); ph = np.unwrap(np.sort(np.angle(Z) % (2 * np.pi), axis=1), axis=0)
    d = R - r
    dphi = ph - ph.mean(0); coll = dphi.mean(1)
    st = (R * (dphi - coll[:, None])).std(); sn = (d - d.mean(0)).std()
    return d.mean(), st, sn
L = []
for k in (2, 4, 6, 8, 10):
    N = 2 * k; lab = f"circle_R31.83_k{k:02d}"; m = M[lab]
    L.append(f"k={k:2d} N={N:2d} measured: depth {np.mean(m['depth']):.2f}, tangential internal {np.sqrt(np.mean(np.square(m['sig_t_internal']))):.2f}, "
             f"normal {np.sqrt(np.mean(np.square(m['sig_n']))):.2f}")
    print(L[-1])
    for K in (32.5, 23.0, 17.8):
        for off in (0.0, 0.61):
            Z, acc = sample(N, K, R - off, sweeps=20000 if N > 12 else 30000)
            for rd in (False, True):
                dm, st, sn = measure(Z, N, rd)
                L.append(f"      K {K:4.1f}, wall offset {off:.2f}, {'grid-rounded' if rd else 'exact       '}: depth {dm:.2f}, tangential {st:.2f}, normal {sn:.2f}  (acc {acc:.2f})")
                print(L[-1])
open(os.path.join(H, "mc_exact.txt"), "w").write("\n".join(L) + "\n")
