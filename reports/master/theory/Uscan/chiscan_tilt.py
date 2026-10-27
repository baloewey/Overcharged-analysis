#!/usr/bin/env python3
"""chi-scan tilt test (2026-10-27). Circles R 31.83 (1.43 l_a), k 4 / 8 / 12, shearSusc chi 0.25 / 0.5 / 1, act 0.07, 3 seeds (runs/chiscan).
Measured exactly as tiltcrowd/depth_correction.one (defect_flow: depth d, coherent-layer tilt; wallshear: wall shear s = e_band + w_band).
Alignment length l_E = sqrt(K / (gamma1 chi s sqrt(lambda^2 - 1))), K 44, gamma1 60, lambda 2.
Question: does the crowding drop of the tilt follow d / l_E (i.e. d sqrt(chi s)), so that chi moves runs along one line?
-> chiscan_tilt.csv, chiscan_tilt.txt"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, "..", "tiltcrowd"))
import depth_correction as DC  # noqa
ROOT = "/home/baloewe/Physics/Overcharged/runs/chiscan"; K, G1, LAM = 44.0, 60.0, 2.0


def one(run):
    r = DC.one(run)
    if not r: return None
    chi = float(run.split("_chi")[1].split("_")[0]); lE = np.sqrt(K / (G1 * chi * r["s"] * np.sqrt(LAM ** 2 - 1)))
    return dict(run=os.path.basename(run), k=r["k"], chi=chi, a=r["a"], d=r["d"], s=r["s"], tilt=r["tilt"], lE=lE, dlE=r["d"] / lE)


if __name__ == "__main__":
    runs = sorted(glob.glob(ROOT + "/circ_k*_chi*_s*"))
    with Pool(int(os.environ.get("NPROC", 14))) as p: res = p.map(one, runs)
    bad = [os.path.basename(r) for r, x in zip(runs, res) if not x or not np.isfinite(x["tilt"])]; res = [x for x in res if x and np.isfinite(x["tilt"])]
    with open(os.path.join(H, "chiscan_tilt.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    L = ["CHI SCAN: coherent-layer tilt vs depth / alignment length (circles R 31.83, act 0.07; mean +- s.e. over seeds)", f"failed / no coherent layer: {bad}",
         "k | chi | seeds | spacing a | depth d | wall shear s | l_E | d/l_E | tilt (deg)"]
    for k in (4, 8, 12):
        for chi in (0.25, 0.5, 1.0):
            g = [x for x in res if x["k"] == k and x["chi"] == chi]
            if not g: continue
            f = lambda c: (np.mean([x[c] for x in g]), np.std([x[c] for x in g], ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.nan)
            L.append(f"{k:2d} | {chi:.2f} | {len(g)} | {f('a')[0]:.1f} | {f('d')[0]:.2f} | {f('s')[0]:.3f} | {f('lE')[0]:.2f} | {f('dlE')[0]:.2f} | {f('tilt')[0]:5.1f} +- {f('tilt')[1]:.1f}")
    x = np.array([r["dlE"] for r in res]); y = np.array([r["tilt"] for r in res]); B, t0 = np.polyfit(x, y, 1)
    L.append(f"joint line tilt = {t0:.1f} + {B:.1f} d/l_E over {len(res)} runs; rms residual {np.std(y - (t0 + B * x)):.1f} deg (raw sd {y.std():.1f})")
    for chi in (0.25, 0.5, 1.0):
        m = np.array([r["chi"] == chi for r in res]); r_ = y[m] - (t0 + B * x[m]); L.append(f"   chi {chi}: mean residual {r_.mean():+.1f} deg")
    open(os.path.join(H, "chiscan_tilt.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
