#!/usr/bin/env python3
"""Tilt potential across crowding, robustness and MPCD-like bulk (user, 2026-10-16; follows k_scan.py). R 31.83, k 3-20, rigid cores, cores fixed.
(1) Numerical robustness: h 0.35 / a 1.5 and h 0.5 / a 2.5 for the point bulk and the split bulk (+1/2 pair at +-0.45 R).
(2) MPCD-like bulk: the MPCD interior is disordered (+1 split, dynamic pairs; ordered only in a ring next to the layer), so the bulk is modelled as 3 +1/2 and
    1 -1/2 at random positions in r < 0.6 R (min separation 4 cells), 6 random configurations, harmonics averaged (positions fixed while the ring registry turns).
(3) Drive calibration (master torque balance; torque on the registry sin tau (D + A1 + 4 A2 cos tau) = 0 -> cos tau* = -(D + A1)/(4 A2), stable if A2 > 0):
    the drive is fixed at k = 6 from the measured tilt (63.8 deg) and extrapolated two ways: per defect constant (D = N d) or per wall length constant
    (D = const, since N a = 2 pi R). Compared with the measured coherent-layer tilt at R 31.83 (bulk_spiral.txt).  -> k_scan2.txt, k_scan2.npz"""
import os, sys, time, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
spec_src = open(os.path.join(H, "relax_frank2.py")).read().split('if __name__ == "__main__":')[0]; RF = {"__file__": os.path.join(H, "relax_frank2.py")}; exec(compile(spec_src, "relax_frank2.py", "exec"), RF)
cfg0 = RF["config"]; R0 = 31.830989
def make_cfg(bulk, seed=0):
    def cfg(k, R, tau, b="point", s=3.0):
        zs, qs = cfg0(k, R, tau, "point")
        zs, qs = zs[:-1], qs[:-1]                                    # wall ring only
        if bulk == "point": zs.append(0j); qs.append(1.0)
        elif bulk == "pair45": zs += [0.45 * R + 0j, -0.45 * R + 0j]; qs += [0.5, 0.5]
        else:
            rng = np.random.default_rng(seed); pts = []
            while len(pts) < 4:
                r, a = 0.6 * R * np.sqrt(rng.random()), 2 * np.pi * rng.random(); z = r * np.exp(1j * a)
                if all(abs(z - p) > 4 for p in pts): pts.append(z)
            zs += pts; qs += [0.5, 0.5, 0.5, -0.5]
        return zs, qs
    return cfg
def one(job):
    k, var, h, a, seed = job; RF["config"] = make_cfg(var, seed)
    r = RF["run"](k=k, R=R0, a=a, h=h, ntau=16); A = RF["harm"](r[:, 2]); return dict(k=k, var=var, h=h, a=a, seed=seed, A=A, conv=int(r[:, 5].sum()))
KS = (3, 4, 5, 6, 8, 10, 12, 15, 20)
MEAS = {}
for l in open(os.path.join(H, "..", "tiltcrowd", "bulk_spiral.txt")):
    p = [c.strip() for c in l.split("|")]
    if len(p) > 4 and p[0].isdigit(): MEAS[int(p[0])] = float(p[3])
if __name__ == "__main__":
    J = [(k, v, h, a, 0) for v in ("point", "pair45") for (h, a) in ((0.35, 1.5), (0.5, 2.5)) for k in KS] + [(k, "rand", 0.5, 1.5, s) for s in range(6) for k in KS]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, J, chunksize=1)
    np.savez(os.path.join(H, "k_scan2.npz"), res=np.array(res, dtype=object), allow_pickle=True)
    groups = {}
    for r in res: groups.setdefault((r["var"], r["h"], r["a"]), {}).setdefault(r["k"], []).append(r["A"])
    L = ["TILT POTENTIAL ACROSS CROWDING, robustness + MPCD-like bulk (R 31.83; whole-ring harmonics A1, A2 in kBT; rand = mean (sd) over 6 random interiors)"]
    for key in sorted(groups):
        L.append(f"-- bulk {key[0]}, h {key[1]}, core {key[2]}"); G = groups[key]
        for k in KS:
            A = np.array(G[k]); m = A.mean(0); sd = A.std(0) if len(A) > 1 else np.zeros(3)
            L.append(f"   k {k:2d} (spacing {2*np.pi*R0/(2*k-2):4.1f}): A1 {m[0]:+8.2f} ({sd[0]:5.2f})  A2 {m[1]:+7.2f} ({sd[1]:5.2f})  per defect A1 {m[0]/(2*k-2):+6.2f}  A2 {m[1]/(2*k-2):+6.3f}")
        # drive calibration and prediction
        A6 = np.array(G[6]).mean(0); tm = np.radians(MEAS.get(6, 63.8)); N6 = 10
        D6 = -A6[0] - 4 * A6[1] * np.cos(tm)
        for lab, Dk in (("drive per defect constant", lambda k: D6 * (2 * k - 2) / N6), ("drive per wall length constant", lambda k: D6)):
            pred = []
            for k in KS:
                A = np.array(G[k]).mean(0); c = -(Dk(k) + A[0]) / (4 * A[1]) if A[1] > 0 else np.nan
                pred.append("slides" if (np.isfinite(c) and c < -1) else ("untilted" if (np.isfinite(c) and c > 1) else (f"{np.degrees(np.arccos(c)):5.1f}" if np.isfinite(c) else "A2<=0")))
            L.append(f"   prediction ({lab}; D6 = {D6:+.1f}): " + "  ".join(f"k{k}:{p}" for k, p in zip(KS, pred)))
        L.append("   measured: " + "  ".join(f"k{k}:{MEAS.get(k, np.nan):5.1f}" for k in KS))
    print("\n".join(L)); open(os.path.join(H, "k_scan2.txt"), "w").write("\n".join(L) + "\n")
