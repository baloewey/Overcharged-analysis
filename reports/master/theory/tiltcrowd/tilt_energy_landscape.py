#!/usr/bin/env python3
"""Elastic energy of tilting crowded wall layers (user, 2026-11-23). Circle R 31.83 (MPCD k-scan size), wall ring N = 2k-2 +1/2 at the depth-law radius
plus a central +1 (crowd_elastic.ring), one Frank constant K = 44 kBT (the revised static constant), rigid cores radius 1.5, field relaxed
(ejection/tilted_structures.energy; a core-director rotation by delta/2 turns a +1/2 comet by delta). Two ways to realise a tilt:
  REGISTRY: positions slide along the wall (rotation of the ring by tau/N), cores at their relaxed (lock) orientation; energy = exact one-constant
            point-defect energy with wall images (pi K analytic_E; the rigid-core grid energy has a few-kBT rasterisation noise when cores move);
  TWIST:    positions fixed, cores rotated to the target tilt (core offset (tau - natural tilt)/2).
(1) COLLECTIVE: all wall defects at a common tilt tau (-90..90 deg): E_reg(tau) (zero mode for one constant up to the central +1 and discretisation), E_tw(tau).
(2) MEAN FIELD: all wall defects at a common tilt tau0 (registry realisation; tau0 = 0 and 60 deg) except defect 0, whose tilt tau1 is varied
    (a) by sliding it alone along the wall (its actual tilt taken from natural_tilts) or (b) by rotating its core at fixed position.
k 4, 6 (sparse controls), 12, 16, 20 (spacing 8.7 -> 5.1 cells). Energies relative to the reference configuration, kBT; collective also per defect.
-> tilt_energy_landscape.txt, tilt_energy_landscape.npz"""
import os, sys, numpy as np
from multiprocessing import Pool
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path.insert(0, T_ + "/ejection")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    exec(compile(open(path).read().split('\nif __name__ == "__main__":')[0], path, "exec"), m.__dict__); return m
TS = load("ts", T_ + "/ejection/tilted_structures.py")
H = os.path.dirname(os.path.abspath(__file__)); R = 31.830989; K = 44.0; KS = (4, 6, 12, 16, 20)
TAUS = np.radians(np.arange(-90, 91, 15)); DT1 = np.radians(np.arange(-90, 91, 15))


def ring(k, al=0.0):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [0j], [0.5] * N + [1.0], N, rs


def E(zs, qs, c): return TS.energy(zs, qs, c, R, K, K)                 # grid, rigid cores: used only at FIXED positions (twist), where rasterisation cancels
def EA(zs, qs): return np.pi * K * TS.analytic_E(zs, qs, R)                 # exact one-constant point-defect energy with wall images: used when positions move


def job(j):
    kind, k, p = j; N = 2 * k - 2
    if kind == "coll_reg":
        zs, qs, _, _ = ring(k, p / N); return (j, EA(zs, qs), float(np.degrees(np.mean(np.angle(np.exp(1j * TS.natural_tilts(zs, qs, R, N)))))))
    if kind == "coll_tw":
        zs, qs, _, _ = ring(k, 0.0); t0 = TS.natural_tilts(zs, qs, R, N); c = list(np.angle(np.exp(1j * (p - t0))) / 2) + [0.0]; return (j, E(zs, qs, c), float(np.degrees(p)))
    tau0, dt = p
    zs, qs, _, rs = ring(k, tau0 / N)
    if kind == "mf_reg":                                                      # slide defect 0 alone by dt/N (its lock tilt moves by ~dt)
        zs = list(zs); zs[0] = rs * np.exp(1j * (dt / N)) * np.exp(1j * (tau0 / N)); t1 = TS.natural_tilts(zs, qs, R, N)[0]
        return (j, EA(zs, qs), float(np.degrees(np.angle(np.exp(1j * t1)))))
    if kind == "mf_tw":                                                       # rotate defect 0's core by dt at fixed positions
        t0 = TS.natural_tilts(zs, qs, R, N); c = [0.0] * (N + 1); c[0] = dt / 2; return (j, E(zs, qs, c), float(np.degrees(np.angle(np.exp(1j * (t0[0] + dt))))))


if __name__ == "__main__":
    jobs = [("coll_reg", k, t) for k in KS for t in TAUS] + [("coll_tw", k, t) for k in KS for t in TAUS]
    jobs += [(m, k, (np.radians(t0), dt)) for m in ("mf_reg", "mf_tw") for k in KS for t0 in (0, 60) for dt in DT1]
    with Pool(int(os.environ.get("NPROC", 16))) as pool: res = pool.map(job, jobs, chunksize=1)
    D = {}
    for (kind, k, p), e, t in res: D.setdefault((kind, k, p if np.isscalar(p) else round(float(np.degrees(p[0])))), []).append((t, e))
    L = ["ELASTIC COST OF TILTING THE WALL LAYER (circle R 31.83, ring 2k-2 + central +1, one constant K 44 kBT, rigid cores 1.5; energies in kBT)"]
    for k in KS:
        a = 2 * np.pi * ring(k)[3] / (2 * k - 2); N = 2 * k - 2; L.append(f"\nk {k} (N {N}, spacing {a:.1f} cells)")
        for kind, lab in (("coll_reg", "COLLECTIVE registry (ring slides): tilt -> E - E(0), per defect"), ("coll_tw", "COLLECTIVE twist (cores rotate in place): tilt -> E - E(0), per defect")):
            rows = sorted((t, e) for p in TAUS for (t, e) in D[(kind, k, p)])
            ez = min(rows, key=lambda x: abs(x[0]))[1]
            L.append(f"  {lab}: " + "  ".join(f"{t:+.0f}:{(e - ez) / N:+.2f}" for t, e in rows))
        for kind, lab in (("mf_reg", "MEAN FIELD, odd defect slides alone"), ("mf_tw", "MEAN FIELD, odd defect rotates in place")):
            for t0 in (0, 60):
                rows = sorted(D[(kind, k, t0)]); ez = min(rows, key=lambda x: abs(x[0] - t0))[1]
                L.append(f"  {lab}, others at {t0} deg: tau1 -> E - E(tau1 = {t0}) " + "  ".join(f"{t:+.0f}:{e - ez:+.1f}" for t, e in rows))
    np.savez(os.path.join(H, "tilt_energy_landscape.npz"), res=np.array([(r[0][0], r[0][1], float(np.degrees(r[0][2])) if np.isscalar(r[0][2]) else float(np.degrees(r[0][2][0])),
                                                                      0.0 if np.isscalar(r[0][2]) else float(np.degrees(r[0][2][1])), r[1], r[2]) for r in res], dtype=object))
    open(os.path.join(H, "tilt_energy_landscape.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
