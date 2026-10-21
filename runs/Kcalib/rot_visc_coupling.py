#!/usr/bin/env python3
"""Rotational viscosity vs flow coupling (2026-10-21): rot_visc.py's mode-decay measurement applied to runs/Kcalib_lambda (shearSusc 0; tumble 0, 1, -2)
and runs/Kcalib L64 (production). Decay rate Gamma_q of splay (q || y) and bend (q || x) modes, m = 2, 3, 4, 6 (m 1 too slow for 1800 t.u.);
gamma_eff = K q^2 / Gamma with each run's own fitted K1/K3 (director-frame ML, reports/master/data/frank_cache.pkl). The decay rate is set by the
deterministic relaxation (response K / friction), not by the noise amplitude, so Gamma_splay/Gamma_bend also probes the response anisotropy
(confounded by backflow when the coupling is on). -> rot_visc_coupling.txt"""
import pickle, numpy as np
from multiprocessing import Pool
from rot_visc import thetas, acf
C = pickle.load(open("/home/baloewe/Physics/Overcharged/reports/master/data/frank_cache.pkl", "rb")); R = "/home/baloewe/Physics/Overcharged/runs/"
RUNS = ["Kcalib/L64_s1", "Kcalib/L64_s2"] + [f"Kcalib_lambda/L64_{c}_s{s}" for c in ("susc0", "tum0", "tum1", "tumm2") for s in (1, 2)]
def one(r):
    TQ, L = thetas(R + r); K1, K3 = C[r]["K1"], C[r]["K3"]; rows = []
    for m in (2, 3, 4, 6):
        q = 2 * np.pi * m / L; G = {}
        for lab, idxs in (("splay", [(0, m), (0, -m)]), ("bend", [(m, 0), (-m, 0)])):
            c = np.mean([acf(TQ[:, i % L, j % L], 60) for i, j in idxs], axis=0); t = np.arange(len(c)) * 2.0
            sel = (c > 0.2) & (t > 0); G[lab] = -np.polyfit(t[sel], np.log(c[sel]), 1)[0]
        rows.append((m, q, G["splay"], G["bend"], K1 * q * q / G["splay"], K3 * q * q / G["bend"]))
    return r, K1, K3, rows
if __name__ == "__main__":
    with Pool(5) as p: res = p.map(one, RUNS)
    L = ["run | K1 K3 (fluct) | per m: Gamma_splay Gamma_bend ratio | gamma_splay gamma_bend | means over m 2-6"]
    for r, K1, K3, rows in res:
        a = np.array(rows); L.append(f"{r}: K1 {K1:.1f} K3 {K3:.1f}")
        for m, q, gs, gb, ys, yb in rows: L.append(f"   m {m} q {q:.3f}: Gamma splay {gs:.4f} bend {gb:.4f} ratio {gs/gb:.2f} | gamma splay {ys:5.1f} bend {yb:5.1f}")
        L.append(f"   MEAN: Gamma ratio splay/bend {np.mean(a[:,2]/a[:,3]):.2f}; gamma splay {a[:,4].mean():.1f}, bend {a[:,5].mean():.1f}")
    open("rot_visc_coupling.txt", "w").write("\n".join(L) + "\n"); print("\n".join(l for l in L if "MEAN" in l or ":" in l and "m " not in l))
