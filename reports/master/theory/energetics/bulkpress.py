#!/usr/bin/env python3
"""Elastic pressure of the bulk charge on wall defects, tips vs flanks (user, 2026-10-12: if the excess charge sits near the middle, it presses less on the tip wall defects).
Shape sweep, second half, every 5th frame. For each wall +1/2 (< 4 cells): push from the bulk defects (> 4 cells) with the 2D defect interaction, f = sum_j 2 q_i q_j (r_i - r_j)/|r_i - r_j|^2
(units of pi K; repulsive for like charges), projected on the OUTWARD normal (> 0 pushes the wall defect towards the wall), and its depth. By aspect and region
(tip |cos t| > 0.866, flank < 0.5, mid otherwise; ellipse angle parameter). No image charges (a proxy). -> bulkpress.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from tipcount import reg  # noqa
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; ainv = bc.get("aInv", [1, 1, 1]); B = float(bc["R"]); A = B * float(ainv[0] or 1); k = bc["customAnch"]; sg = 1 if k > 0 else -1
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5); out = []
        for f in fr[len(fr) // 2::5]:
            if not len(f): continue
            d, j = tree.query(f[:, :2]); w = np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)); b = np.flatnonzero(d > 4)
            for i in w:
                zi = f[i, 0] + 1j * f[i, 1]; push = 0.0
                for jj in b:
                    dz = zi - (f[jj, 0] + 1j * f[jj, 1]); push += 2 * 0.5 * np.sign(f[i, 2]) * 0.5 * np.sign(f[jj, 2]) * (dz / abs(dz) ** 2 * np.conj(Nn[j[i]])).real
                out.append((reg(f[i, :2], A, B), push, d[i], 0.5 * np.sign(f[b, 2]).sum() if len(b) else 0.0))
        return (A / B, np.array(out).reshape(-1, 4))
    except Exception as ex:
        return (None, run + " " + str(ex))
if __name__ == "__main__":
    SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; runs = sorted(glob.glob(SW + "/ellipse/*/*/*/s0*") + glob.glob(SW + "/circle/circle/*/*/s0*"))
    with Pool(int(os.environ.get("NPROC", 3))) as pool: res = [r for r in pool.map(one, runs, chunksize=2) if r[0] is not None]
    L = ["Push of the bulk charge on wall +1/2 along the outward normal (> 0 = towards the wall; units pi K per cell), and depth, by aspect and region. Mean ± s.e. over samples.",
         "aspect | region | n | push | depth | push per unit bulk charge (slope within aspect)"]
    for a in sorted(set(round(r[0], 2) for r in res)):
        X = np.concatenate([r[1] for r in res if round(r[0], 2) == a and len(r[1])])
        for g, nm in ((0, "tip"), (1, "mid"), (2, "flank")):
            m = X[:, 0] == g
            if m.sum() < 50: continue
            sl = np.polyfit(X[m, 3], X[m, 1], 1)[0] if np.std(X[m, 3]) > 0 else np.nan
            L.append(f"  {a:4.2f} | {nm:5s} | {m.sum():6d} | {X[m, 1].mean():+.4f} ± {X[m, 1].std() / np.sqrt(m.sum()):.4f} | {X[m, 2].mean():.2f} | {sl:+.4f}")
    print("\n".join(L)); open(os.path.join(H, "bulkpress.txt"), "w").write("\n".join(L) + "\n")
