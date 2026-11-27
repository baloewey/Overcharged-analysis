#!/usr/bin/env python3
"""Lock slope from fluctuations in DRIVEN PASSIVE cavities (user, 2026-11-27): same estimator as lock_slope_fluct.py (frames with exactly 2k-2 wall +1/2,
dtau_j vs dphi_j about the frame's regular ring, OLS slope vs 2k-2), on imposed-vortex runs: swirl_thin (Test 2: thin / medium / disk forcing),
swirl_crowd A03 (disk forcing, k 6-16), swirl_shear (k 6 / 13, swirl amplitude scan), swirl_threshold (k 6, amplitude scan).
Also the common tilt (circular mean of tau in those frames) and the displacement spread. -> lock_slope_driven.txt"""
import glob, os, sys, re, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import lock_slope_fluct as LF  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs"
SETS = [("Test2 thin", R_ + "/swirl_thin/thin_k*_s*"), ("Test2 medium", R_ + "/swirl_thin/medium_k*_s*"), ("Test2 disk", R_ + "/swirl_thin/disk_k*_s*"),
        ("Disk A03", R_ + "/swirl_crowd/A03_k*_s*"), ("Shear scan", R_ + "/swirl_shear/k*_A*"), ("Threshold k6", R_ + "/swirl_threshold/passive_A*_s*")]


def tag(run):
    for lab, pat in SETS:
        if run in glob.glob(pat): return lab
    return "?"


if __name__ == "__main__":
    runs = [d for _, pat in SETS for d in sorted(glob.glob(pat)) if os.path.isdir(d) and (glob.glob(d + "/defects.dat*"))]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(LF.run_one, runs, chunksize=1)
    rows = {}
    for run, r in zip(runs, res):
        if not r or "err" in r: continue
        lab = tag(run); extra = ""
        if lab in ("Shear scan", "Threshold k6"): m = re.search(r"A([0-9.]+)", os.path.basename(run)); extra = f" A{m.group(1).rstrip('.')}" if m else ""
        rows.setdefault((lab + extra, r["k"]), []).append(r)
    L = ["LOCK SLOPE FROM FLUCTUATIONS, DRIVEN PASSIVE CAVITIES (frames with exactly 2k-2 wall +1/2; second half)",
         "set | k | defect-frames | sd displacement (cells) | OLS slope | 2k-2 | slope/(2k-2) | corr"]
    for (lab, k) in sorted(rows, key=lambda x: (x[0], x[1])):
        g = rows[(lab, k)]; x = np.concatenate([r["x"] for r in g]); y = np.concatenate([r["y"] for r in g])
        if len(x) < 200: continue
        sl = np.cov(x, y)[0, 1] / np.var(x); cr = np.corrcoef(x, y)[0, 1]; R0 = 31.83
        L.append(f"{lab:18s} | {k:2d} | {len(x):6d} | {x.std() * R0:4.2f} | {sl:6.2f} | {2 * k - 2:2d} | {sl / (2 * k - 2):.2f} | {cr:.2f}")
    bad = [run for run, r in zip(runs, res) if r and "err" in r][:3]
    for b in bad: L.append("  ERROR " + b)
    open(os.path.join(H, "lock_slope_driven.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
