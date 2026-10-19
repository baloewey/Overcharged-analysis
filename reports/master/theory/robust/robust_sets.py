#!/usr/bin/env python3
"""Paper-relevant run sets (analysis item 2, 2026-10-19), measured with energetics/tiltlaw.one (second half: wall +1/2 < 4 cells, bulk charge >= 4 cells,
wall-defect tilt |tau| median/mean, tilt coherence |<e^(i tau)>| per frame, |wall v_t|, v_rms, P(circulation keeps sign)), grouped by run name without the seed.
Sets: act_robust (zeta 0.06/0.08, circle k 6 R/l_a 1.43 and ellipse eps 2 R_eq 1.5 l_a k 7, sizes scaled by l_a ~ zeta^-1/2), lambda_crowd (circle area 2500 k 14,
lambda 2-6) with lambda_circle k 6 (same circle) as the sparse reference, lambda_ellipse (eps 2, area 2500, k 5, lambda 2-6).
Also coherent-layer tilt (frames with C > 0.9, median |arg <e^(i tau)>|) computed here. -> robust_sets.txt, robust_sets.pkl"""
import os, sys, glob, json, re, pickle, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/energetics")
import tiltlaw as TL  # noqa
import energetics as En  # noqa
from scipy.spatial import cKDTree
RUNS = "/home/baloewe/Physics/Overcharged/runs"
def coh_tilt(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; k = bc["customAnch"]
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = En.read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); fr = fr[len(fr) // 2:]; tc = []; C = []
    for f in fr:
        d, j = tree.query(f[:, :2]) if len(f) else (np.zeros(0), np.zeros(0, int)); w = (f[:, 2] > 0) & (d < 4) if len(f) else np.zeros(0, bool)
        if w.sum() < 3: continue
        z = np.mean(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))); C.append(abs(z))
        if abs(z) > 0.9: tc.append(abs(np.angle(z)))
    return np.degrees(np.median(tc)) if len(tc) > 20 else np.nan, np.mean(np.array(C) > 0.9) if C else np.nan
def one(run):
    try:
        r = TL.one(run); ct, fc = coh_tilt(run); return run, r, ct, fc
    except Exception as e:
        return run, None, repr(e), None
if __name__ == "__main__":
    sets = {"act_robust": "act_robust/*_s?", "lambda_crowd": "lambda_crowd/k14_*", "lambda_circle_k6": "lambda_circle/k06_*", "lambda_ellipse": "lambda_ellipse/eps2_*"}
    runs = [(nm, r) for nm, g in sets.items() for r in sorted(glob.glob(f"{RUNS}/{g}")) if os.path.exists(r + "/DONE")]
    with Pool(6) as p: res = p.map(one, [r for _, r in runs])
    pickle.dump(res, open("robust_sets.pkl", "wb"))
    G = {}
    for (nm, run), (_, r, ct, fc) in zip(runs, res):
        if r is None: print("FAILED", run, ct); continue
        G.setdefault((nm, re.sub(r"_s\d+$", "", os.path.basename(run))), []).append(list(r[2:11]) + [ct, fc] + [r[0]])
    L = ["set | group | n | k | wall +1/2 (2k-2) | bulk charge | tilt |tau| median / mean | arccos(1/lambda) | coherence C | coherent-layer tilt (C>0.9) | frac C>0.9 | |wall v_t| | v_rms | P(sign kept)"]
    for key in sorted(G):
        S = np.array(G[key], float); m = np.nanmean(S, 0); e = np.nanstd(S, 0) / np.sqrt(max(len(S) - 1, 1)); k = int(S[0, -1])
        lam = re.search(r"lam([\d.]+)", key[1]); lam = float(lam.group(1)) if lam else 2.0
        L.append(f"{key[0]:16s} | {key[1]:20s} | {len(S)} | {k:2d} | {m[0]:5.2f}±{e[0]:.2f} ({2*k-2}) | {m[1]:+.2f}±{e[1]:.2f} | {m[2]:5.1f} / {m[3]:5.1f}±{e[3]:.1f} | {np.degrees(np.arccos(1/lam)):5.1f} | "
                 f"{m[4]:.2f} | {m[9]:5.1f}±{e[9]:.1f} | {m[10]:.2f} | {m[5]:.3f} | {m[6]:.3f} | {m[7]:.2f}")
    open("robust_sets.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
