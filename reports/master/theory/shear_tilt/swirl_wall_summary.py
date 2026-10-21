#!/usr/bin/env python3
"""Shear-only control with the force over the whole disk (tag _wall, 2026-10-20): laminar zeta 3 disk + azimuthal body force A everywhere, k 4 vs 16.
Same estimators as zeta_scan_summary.py (wall shear s over 1-4 cells, d/l_E). -> swirl_wall_summary.txt"""
import json, os, numpy as np
import shear_tilt as M
K, G1 = 1.535, 1.53; L = ["SHEAR-ONLY CONTROL, force over the whole disk: laminar zeta 3 + azimuthal body force A; target zeta 8: wall shear 1.87 (k 4) / 2.14 (k 16), drop 10.7.",
    "A k | tilt | coh | circ | nwall (2k-2) | bulk +/- | spiral | tilt sd | d | wall shear s | d/l_E | drop k4-k16"]
T = {}
for A in (0.02, 0.05, 0.1, 0.2):
    for k in (4, 16):
        fn = f"var_swirl_{A:g}_k{k:02d}_wall.json"
        if not os.path.exists(fn): continue
        d = json.load(open(fn)); ts = [t for t in d.get("tilt_series", []) if t is not None]; sd = np.std(ts[len(ts) // 2:]) if ts else np.nan; T[(A, k)] = d["tilt"]
        dle = "| - | - | -"; ff = f"fields_swirl_{A:g}_k{k:02d}_wall.npz"
        if os.path.exists(ff):
            F = np.load(ff); q1, q2, vx, vy = F["q1"], F["q2"], F["vx"], F["vy"]; R = float(F["R"]); M.k = k
            D = M.defects(q1, q2); dd = np.array([R - np.hypot(px, py) for px, py, q, a in D if q > 0 and np.hypot(px, py) > R - 6])
            r = np.hypot(M.X, M.Y); ph = np.arctan2(M.Y, M.X); vt = -vx * np.sin(ph) + vy * np.cos(ph); e = np.arange(R - 8, R + 0.01, 0.5); rc = 0.5 * (e[1:] + e[:-1])
            prof = np.array([vt[M.inside & (r >= a) & (r < b)].mean() for a, b in zip(e[:-1], e[1:])]); b = (R - rc >= 1) & (R - rc <= 4)
            s = abs(np.polyfit(rc[b], prof[b], 1)[0] - np.mean(prof[b] / rc[b])); lE = np.sqrt(K / (G1 * s * np.sqrt(3)))
            if len(dd): dle = f"| {np.median(dd):.2f} | {s:.3f} | {np.median(dd)/lE:.2f}"
            else: dle = f"| - | {s:.3f} | -"
        tl = d["tilt"]; drop = (f"| {T[(A,4)]-T[(A,16)]:.1f}" if k == 16 and T.get((A, 4)) is not None and tl is not None and d.get('coh', 0) > 0.9 else "|")
        L.append(f"{A:<5g} {k:2d} | {'-' if tl is None else f'{tl:.1f}'} | {d.get('coh', np.nan):.2f} | {d['circ']:+.2f} | {d['nwall']:.1f} ({2*k-2}) | {d['nbulk_p']:.1f}/{d['nbulk_m']:.1f} | {d['spiral_order']:.2f} | {sd:.1f} {dle} {drop}")
open("swirl_wall_summary.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
