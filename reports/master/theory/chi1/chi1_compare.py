#!/usr/bin/env python3
"""chi = 1 vs chi = 0.5 MPCD circles at matched R/l_a and k (user, 2026-11-23/24). chi 1: runs/chi1_matched (l_a 15.1); chi 0.5: main sweep
runs/cavity/server_main/main_act0.07 (l_a 22.2), runs with R/l_a within 0.08 of the chi-1 value and the same k. Same estimator for both (defects.dat,
second half): wall +1/2 (< 4 cells), bulk charge, layer tilt |arg<exp(i tau)>| over frames with tilt coherence C > 0.9 (as the tilt-law tables),
fraction of coherent frames. Spacing a = 2 pi r_layer/(2k-2) ~ 2 pi R/(2k-2) in cells and in units of l_a.
-> chi1_compare.txt"""
import glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); RUNS = "/home/baloewe/Physics/Overcharged/runs"


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
    nw, qb, tl, nco = [], [], [], 0
    for f in fr:
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4); nw.append(w.sum()); qb.append(0.5 * np.sign(f[R - r >= 4, 2]).sum())
        if w.sum() >= 3:
            z = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))
            if abs(z) > 0.9: tl.append(z); nco += 1
    tilt = abs(np.degrees(np.angle(np.mean(tl)))) if tl else np.nan
    return dict(run=run, R=R, k=k, nw=np.mean(nw), qb=np.mean(qb), tilt=tilt, fco=nco / max(len(fr), 1))


if __name__ == "__main__":
    c1 = sorted(os.path.dirname(p) for p in glob.glob(RUNS + "/chi1_matched/R*_k*_s*/DONE"))
    targets = {}
    for r in c1:
        rl = float(os.path.basename(r).split("_")[0][1:]); k = int(os.path.basename(r).split("_k")[1][:2]); targets[(rl, k)] = True
    c05 = []
    for p in glob.glob(RUNS + "/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/defects.dat.gz"):
        d = os.path.dirname(p); bc = json.load(open(d + "/input.json"))["BC"][0]; rl = bc["R"] / 22.21
        for (tr, tk) in targets:
            if bc["customAnch"] == tk and abs(rl - tr) < 0.08: c05.append((d, tr)); break
    with Pool(12) as p: A = p.map(one, c1); B = p.map(one, [d for d, _ in c05])
    for b, (_, tr) in zip(B, c05): b["rl"] = tr
    for a in A: a["rl"] = float(os.path.basename(a["run"]).split("_")[0][1:])
    L = ["CHI 1 vs CHI 0.5 at matched R/l_a and k (MPCD circles, lambda 2, act 0.07). Means over seeds.",
         "R/l_a | k | spacing (cells chi1 / chi0.5; units of l_a) | chi 1: wall (2k-2), bulk, tilt (coherent frac) | chi 0.5: same | runs"]
    for rl in (1.0, 1.43, 2.0):
        for k in sorted({a["k"] for a in A if a["rl"] == rl}):
            ga = [a for a in A if a["rl"] == rl and a["k"] == k]; gb = [b for b in B if b["rl"] == rl and b["k"] == k]
            f = lambda g, key: np.nanmean([x[key] for x in g]) if g else np.nan
            sa = 2 * np.pi * f(ga, "R") / (2 * k - 2); sb = 2 * np.pi * f(gb, "R") / (2 * k - 2) if gb else np.nan
            L.append(f"{rl:4.2f} | {k:2d} | {sa:5.1f} / {sb:5.1f} ({sa / 15.1:.2f} l_a) | {f(ga, 'nw'):5.2f} ({2 * k - 2}), {f(ga, 'qb'):+.2f}, {f(ga, 'tilt'):5.1f} ({f(ga, 'fco'):.2f}) | "
                     f"{f(gb, 'nw'):5.2f}, {f(gb, 'qb'):+.2f}, {f(gb, 'tilt'):5.1f} ({f(gb, 'fco'):.2f}) | {len(ga)} / {len(gb)}")
    open(os.path.join(H, "chi1_compare.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
