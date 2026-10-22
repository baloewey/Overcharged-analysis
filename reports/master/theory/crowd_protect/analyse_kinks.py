#!/usr/bin/env python3
"""Summary of kink_life.npz by wall density level and shape class. -> analyse_kinks.txt"""
import csv, numpy as np
Z = np.load("kink_life.npz", allow_pickle=True); SR = {r["run"]: r for r in csv.DictReader(open("../shapes/shape_runs.csv"))}
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); runs = Z["runs"]; HB = Z["HB"]; hc = 0.5 * (HB[:-1] + HB[1:])
rho = np.array([float(SR[r]["rho"]) for r in runs]); asp = np.array([float(SR[r]["aspect"]) for r in runs]); lev = np.argmin(np.abs(rho[:, None] - LV), 1)
L = ["KINKS (wrong tilt sign relative to the frame majority), tracked every 2 t.u.; medians over runs unless stated",
     "class | rho | wrong fraction | C | sign changes/ring | episode rate /defect/1000 t.u. | episode duration median / mean / P(>10 t.u.) / P(>40 t.u.) | tilt hist: majority peak, share near -peak (true flip), share beyond +-110 (head-flip)"]
for nm, sel in (("circles", asp < 1.05), ("aspect>=1.4", asp >= 1.4)):
    for i, r0 in enumerate(LV):
        m = np.flatnonzero((lev == i) & sel); F = [Z["fstat"][j] for j in m]
        wf = np.median([1 - f[:, 1].mean() for f in F]); C = np.median([f[:, 0].mean() for f in F]); sc = np.median([f[:, 2].mean() for f in F])
        E = np.concatenate([Z["eps"][j] for j in m]); h = Z["hist"][m].sum(0); h = h / h.sum(); pk = hc[hc > 0][np.argmax(h[hc > 0])]
        near = h[(hc < -pk + 25) & (hc > -pk - 25)].sum(); far = h[np.abs(hc) > 110].sum()
        L.append(f"{nm:11s} | {r0:.2f} | {wf:.3f} | {C:.2f} | {sc:.1f} | {np.median(Z['rate'][m]):5.1f} | {np.median(E):4.0f} / {E.mean():5.1f} / {np.mean(E > 10):.2f} / {np.mean(E > 40):.3f} | "
                 f"{pk:+.0f}, {near:.3f}, {far:.3f}")
open("analyse_kinks.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
