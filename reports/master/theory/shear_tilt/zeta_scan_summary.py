#!/usr/bin/env python3
"""Continuum active disk, activity scan zeta 0.3-10 (k 4 vs 16; 2026-10-20): tilt, coherence, circulation, wall count, bulk +-1/2, spiral order, tilt sd (second half)
from var_act_*.json; core depth d, wall shear s and d/l_E (K 1.535, gamma1 1.53, lambda 2) from the saved fields (fields_act_<z>_k<k>.npz; zeta 3: *_fields.npz).
-> zeta_scan_summary.txt"""
import json, os, numpy as np
import shear_tilt as M
K, G1 = 1.535, 1.53; L = ["zeta k | tilt | coh | circ | nwall (2k-2) | bulk +/- | spiral order | tilt sd 2nd half | d median | wall shear | d/l_E | drop k4-k16"]
T = {}
for z in (3, 4, 5, 6, 7, 8, 10):
    for k in (4, 16):
        fn = f"var_act_{z:g}_k{k:02d}.json"
        if not os.path.exists(fn): continue
        d = json.load(open(fn)); ts = [t for t in d["tilt_series"] if t is not None]; sd = np.std(ts[len(ts) // 2:]) if ts else np.nan; T[(z, k)] = d["tilt"]
        dle = "| - | - | -"
        for ff in (f"fields_act_{z:g}_k{k:02d}.npz", f"fields_act_{z:g}_k{k:02d}_fields.npz"):
            if not os.path.exists(ff): continue
            F = np.load(ff); q1, q2, vx, vy = F["q1"], F["q2"], F["vx"], F["vy"]; R = float(F["R"]); M.k = k
            D = M.defects(q1, q2); dd = np.array([R - np.hypot(px, py) for px, py, q, a in D if q > 0 and np.hypot(px, py) > R - 6])
            r = np.hypot(M.X, M.Y); ph = np.arctan2(M.Y, M.X); vt = -vx * np.sin(ph) + vy * np.cos(ph); e = np.arange(R - 8, R + 0.01, 0.5); rc = 0.5 * (e[1:] + e[:-1])
            prof = np.array([vt[M.inside & (r >= a) & (r < b)].mean() for a, b in zip(e[:-1], e[1:])]); b = (R - rc >= 1) & (R - rc <= 4)
            s = abs(np.polyfit(rc[b], prof[b], 1)[0] - np.mean(prof[b] / rc[b])); lE = np.sqrt(K / (G1 * s * np.sqrt(3)))
            if len(dd): dle = f"| {np.median(dd):.2f} | {s:.3f} | {np.median(dd)/lE:.2f}"
            break
        tl = d["tilt"]; drop = (f"| {T[(z,4)]-T[(z,16)]:.1f}" if k == 16 and (z, 4) in T and T[(z, 4)] is not None and tl is not None and d.get('coh', 0) > 0.9 else "|")
        L.append(f"{z:>4g} {k:2d} | {'-' if tl is None else f'{tl:.1f}'} | {d.get('coh', np.nan):.2f} | {d['circ']:+.2f} | {d['nwall']:.1f} ({2*k-2}) | {d['nbulk_p']:.1f}/{d['nbulk_m']:.1f} | {d['spiral_order']:.2f} | {sd:.1f} {dle} {drop}")
open("zeta_scan_summary.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
