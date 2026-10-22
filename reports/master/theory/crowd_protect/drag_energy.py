#!/usr/bin/env python3
"""Elastic cost of dragging one wall +1/2 out of the layer, vs wall density (user, 2026-10-22: 'energy of dragging a defect with a given tilt against the
winding for dense layers; in the bulk it should approach the Gauss theorem; what happens near the charges?').
Circle R = 31.83 (P = 200), k = 4, 6, 9, 12, 14 (rho = 2k/P = 0.04-0.14, the shape-sweep levels), ring of N = 2k-2 wall +1/2 at the depth-law radius plus a central +1
(bulk charge 1; as lockvalley/lock_valley.py). Layer registry alpha = tau_L/N so the relaxed tilts equal tau_L (the lock). Probe = wall defect 0, moved radially
inward by delta; all others fixed in position.
  FIXED : every wall core (probe included) held at tilt tau_L (flow-fixed tilts; rigid cores a = 1.5, field outside relaxed; tilted_structures.energy)
  FREE  : natural core orientations (c = 0; no imposed tilt)
  one constant K = 44 kBT. Energies relative to delta = 0, kBT.
GAUSS (analytic, point charges, one constant): probe energy vs delta with the other N-1 wall defects discrete vs smeared into a uniform ring of the same charge
at the same radius (2D shell theorem: a smeared ring exerts no force inside it). The difference is the discreteness ('near the charges') part.
-> drag_energy.txt, drag_energy.npz"""
import os, sys, numpy as np
from multiprocessing import Pool
import importlib.util
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    src = open(path).read().split('\nif __name__ == "__main__":')[0]; exec(compile(src, path, "exec"), m.__dict__); return m
TS = load("ts", T_ + "/ejection/tilted_structures.py")
R = 31.830989; K = 44.0; KS = [4, 6, 9, 12, 14]; TAUS = [0, 30, 60, 90]; DEL = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12]
def ring(k, tau):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; al = np.radians(tau) / N
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [0j], [0.5] * N + [1.0], N, rs
def config(k, tau, dl):
    zs, qs, N, rs = ring(k, tau); zs[0] = zs[0] * (rs - dl) / rs; return zs, qs, N
def one(job):
    k, tau, dl, mode = job; zs, qs, N = config(k, tau, dl)
    if mode == "fixed":
        t0 = TS.natural_tilts(zs, qs, R, N); c = list(np.angle(np.exp(1j * (np.radians(tau) - t0))) / 2) + [0.0]
    else: c = [0.0] * (N + 1)
    return job, TS.energy(zs, qs, c, R, K, K)
def probe_E(zp, others, qo):                     # analytic one-constant energy terms involving the probe (q 1/2), units of pi K
    E = -2 * 0.5 * sum(q * np.log(abs(zp - z)) for z, q in zip(others, qo)) - 2 * 0.5 * sum(q * np.log(abs(R ** 2 - zp * np.conj(z))) for z, q in zip(others, qo))
    return E - 0.25 * np.log(R ** 2 - abs(zp) ** 2)
if __name__ == "__main__":
    jobs = [(k, tau, dl, mode) for k in KS for tau in TAUS for dl in DEL for mode in ("fixed", "free")]
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = dict(p.map(one, jobs, chunksize=1))
    L = ["DRAG ENERGY E(delta) - E(0) (kBT, one constant K 44, core 1.5): one wall +1/2 pulled radially inward from the layer, circle R 31.83, central +1"]
    out = {}
    for k in KS:
        N = 2 * k - 2; rs = ring(k, 0)[3]; L.append(f"\nk {k} (N {N}, rho {2 * k / (2 * np.pi * R):.3f}, spacing {2 * np.pi * rs / N:.1f} cells, ring radius {rs:.2f}): columns delta = {DEL}")
        for mode in ("fixed", "free"):
            for tau in (TAUS if mode == "fixed" else [0]):
                e = np.array([res[(k, tau, dl, mode)] for dl in DEL]); e -= e[0]; out[f"{mode}_k{k}_t{tau}"] = e
                L.append(f"  {mode:5s} tau_L {tau:3d}: " + " ".join(f"{v:7.1f}" for v in e) + f" | barrier {e.max():.1f}")
        zs, qs, N, rs = ring(k, 0); oth, qo = zs[1:], qs[1:]
        M = 720; sm = list(rs * np.exp(2j * np.pi * (np.arange(M) + 0.5) / M)); smq = [0.5 * (N - 1) / M] * M
        ed = np.array([probe_E(zs[0] * (rs - dl) / rs, oth, qo) for dl in DEL]); es = np.array([probe_E(zs[0] * (rs - dl) / rs, sm + [0j], smq + [1.0]) for dl in DEL])
        ed = np.pi * K * (ed - ed[0]); es = np.pi * K * (es - es[0]); out[f"gauss_disc_k{k}"] = ed; out[f"gauss_smear_k{k}"] = es
        L.append(f"  analytic discrete : " + " ".join(f"{v:7.1f}" for v in ed)); L.append(f"  analytic smeared  : " + " ".join(f"{v:7.1f}" for v in es))
    open("drag_energy.txt", "w").write("\n".join(L) + "\n"); np.savez("drag_energy.npz", DEL=DEL, **out); print("\n".join(L))
