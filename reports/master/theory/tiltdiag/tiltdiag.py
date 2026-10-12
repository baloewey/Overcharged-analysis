#!/usr/bin/env python3
"""Tilt diagnostics (runs/tiltdiag; passive circle k 6, R 31.83, lam 2, swirl 0.01; density 10/20/40 and mfpot 5/10/20).
Per run, production second half:
 (a) wall +1/2 (depth < 4): count, tilt (tail - normal, coherent mean, signed with the circulation sense of the swirl = +), coherence; bulk charge.
 (b) density per cell vs depth (densityfield), S vs depth and vs distance from the wall-defect cores (directorfield, 2 t.u., matched to defects).
 (c) if detailedSP0.dat exists: cell-level S2, S4 from particle orientations (frames every 50 t.u.), pooled lam_eff = lam (1 + <S4>)/(2 <S2>) by depth and core distance,
     compared with the Maier-Saupe S4 at the same S2.
Prediction: arccos(1/lam_eff) in the region that sets the defect orientation (cells within 1-1.5 of the cores) vs the measured tilt.
Usage: tiltdiag.py RUN [RUN ...] -> tiltdiag_<run>.txt"""
import gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt")
from edge_criteria import read_with_angle  # noqa
from mpcd_profile import S4_of  # noqa
def opener(fn): return gzip.open(fn + ".gz", "rt") if not os.path.exists(fn) and os.path.exists(fn + ".gz") else open(fn)
def frames(fn, ncol, t0, t1, cols):
    cur, tc = [], None
    with opener(fn) as f:
        for line in f:
            p = line.split()
            if len(p) < ncol or not p[0].replace(".", "").replace("e+", "").replace("-", "").isdigit(): continue
            try: t = float(p[0])
            except ValueError: continue
            if t != tc:
                if cur and tc is not None and t0 <= tc <= t1: yield tc, np.array(cur)
                cur, tc = [], t
                if t > t1: return
            if t0 <= t <= t1: cur.append([float(p[c]) for c in cols])
    if cur and tc is not None and t0 <= tc <= t1: yield tc, np.array(cur)
def analyse(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Qc = bc["Q"]; R = bc["R"]; lam = J["species"][0].get("tumble", 2.0)
    dens = J["species"][0]["pop"] / (np.pi * R * R); U = J["species"][0].get("mfpot", 10.0); T = J["simSteps"] * 0.1; t0, t1 = T / 2, T
    L = [f"{os.path.basename(run)}: density {dens:.1f}/cell, mfpot {U}, lam {lam}; window {t0:.0f}-{t1:.0f}"]
    # (a) defects
    td, fr = read_with_angle(run, np.array(Qc[:2], float), 0.0); dmap = {}; nw, z, qb = [], [], []
    for t, f in zip(td, fr):
        if not (t0 <= t <= t1) or not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (r > R - 4) & (f[:, 2] > 0); dmap[round(t, 1)] = f[w, :2]
        nw.append(w.sum()); qb.append(0.5 * np.sign(f[r <= R - 6, 2]).sum())
        if w.any(): z.append(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))))
    zz = np.mean(z); L.append(f"  wall +1/2 {np.mean(nw):.2f} (mode {np.bincount(nw).argmax()}), bulk charge {np.mean(qb):+.2f}, tilt {np.degrees(np.angle(zz)):+.1f} deg, coherence {np.mean(np.abs(z)):.2f}")
    # (b) density vs depth
    if os.path.exists(run + "/densityfield.dat") or os.path.exists(run + "/densityfield.dat.gz"):
        acc = {}
        for k, (t, c) in enumerate(frames(run + "/densityfield.dat", 5, t0, t1, (1, 2, 4))):
            if k % 5: continue
            x = c[:, 0] + 0.5 - Qc[0]; y = c[:, 1] + 0.5 - Qc[1]; d = R - np.hypot(x, y); m = d > -0.5
            for dd, n in zip(d[m], c[m, 2]): b = int(np.floor(dd * 2)) / 2; acc.setdefault(b, []).append(n)
        L.append("  particles per cell vs depth: " + " ".join(f"{b:+.1f}:{np.mean(v):.1f}" for b, v in sorted(acc.items()) if b < 6))
    # S vs depth and near cores (director field)
    accd, accc = {}, {rc: [] for rc in (1.0, 1.5, 2.0, 3.0)}
    for t, c in frames(run + "/directorfield.dat", 8, t0, t1, (1, 2, 7)):
        key = round(t, 1)
        if key not in dmap or not len(dmap[key]): continue
        x = c[:, 0] + 0.5 - Qc[0]; y = c[:, 1] + 0.5 - Qc[1]; d = R - np.hypot(x, y); S = c[:, 2]; m = (d > 0) & (S > 0)
        dc, _ = cKDTree(dmap[key]).query(np.column_stack([x, y]))
        for dd, s in zip(d[m & (dc > 4)], S[m & (dc > 4)]): accd.setdefault(int(dd * 2) / 2, []).append(s)
        for rc in accc: accc[rc].extend(S[m & (dc < rc)])
    L.append("  S vs depth (away from cores): " + " ".join(f"{b:.1f}:{np.mean(v):.3f}" for b, v in sorted(accd.items()) if b < 8))
    L.append("  S near wall cores: " + " ".join(f"r<{rc}: {np.mean(v):.3f}" for rc, v in accc.items() if v))
    # (c) particle-level moments
    fn = run + "/detailedSP0.dat"
    if os.path.exists(fn) or os.path.exists(fn + ".gz"):
        n2 = int(np.ceil(2 * R + 4)); S2d, S4d, S2c, S4c = {}, {}, {rc: [] for rc in (1.0, 1.5, 2.0)}, {rc: [] for rc in (1.0, 1.5, 2.0)}; nfr = 0
        for t, P in frames(fn, 11, t0, t1, (1, 2, 8, 9)):
            key = round(t, 1)
            if key not in dmap: continue
            nfr += 1; ix = np.floor(P[:, 0]).astype(int); iy = np.floor(P[:, 1]).astype(int); psi = np.arctan2(P[:, 3], P[:, 2])
            cid = ix * 1000 + iy; u, inv, cnt = np.unique(cid, return_inverse=True, return_counts=True)
            e2 = np.bincount(inv, np.cos(2 * psi)) + 1j * np.bincount(inv, np.sin(2 * psi)); e2 /= cnt
            th = np.angle(e2) / 2; e4 = np.bincount(inv, np.cos(4 * (psi - th[inv]))) / cnt; S2 = np.abs(e2)
            cx = u // 1000 + 0.5 - Qc[0]; cy = u % 1000 + 0.5 - Qc[1]; d = R - np.hypot(cx, cy); m = (d > 0) & (cnt >= 3)
            dc, _ = cKDTree(dmap[key]).query(np.column_stack([cx, cy])) if len(dmap[key]) else (np.full(len(u), 99.0), None)
            for dd, a, b in zip(d[m & (dc > 4)], S2[m & (dc > 4)], e4[m & (dc > 4)]):
                k2 = int(dd * 2) / 2; S2d.setdefault(k2, []).append(a); S4d.setdefault(k2, []).append(b)
            for rc in S2c: S2c[rc].extend(S2[m & (dc < rc)]); S4c[rc].extend(e4[m & (dc < rc)])
        L.append(f"  particle moments ({nfr} frames): depth: S2, S4 (Maier-Saupe S4 at that S2), lam_eff")
        for k2 in sorted(S2d):
            if k2 >= 6: continue
            a, b = np.mean(S2d[k2]), np.mean(S4d[k2]); L.append(f"     {k2:.1f}: S2 {a:.3f} S4 {b:.3f} (MS {S4_of(min(a, .995)):.3f}) lam_eff {lam * (1 + b) / (2 * a):.2f}")
        for rc in S2c:
            if not S2c[rc]: continue
            a, b = np.mean(S2c[rc]), np.mean(S4c[rc]); le = lam * (1 + b) / (2 * a); ms = lam * (1 + S4_of(min(a, .995))) / (2 * a)
            L.append(f"   cores r<{rc}: S2 {a:.3f} S4 {b:.3f} -> lam_eff {le:.2f} (Maier-Saupe {ms:.2f}) -> arccos(1/lam_eff) {np.degrees(np.arccos(1 / le)):.1f} (MS {np.degrees(np.arccos(1 / ms)):.1f})")
    return L
if __name__ == "__main__":
    for run in sys.argv[1:]:
        L = analyse(run.rstrip("/")); print("\n".join(L), flush=True)
        open(f"tiltdiag_{os.path.basename(run.rstrip('/'))}.txt", "w").write("\n".join(L) + "\n")
