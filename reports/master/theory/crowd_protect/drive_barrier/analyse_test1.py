#!/usr/bin/env python3
"""Analysis of test 1. Near-layer active temperature T_near per run: least-squares fit of ln S(x) = a - U(x)/T_near over x in [d0 + 0.5, d0 + 2.5]
(S from drive_barrier.csv; U = relaxed-neighbour drag energy, drag_relaxed.npz). Drive measures:
  (a) FLOW: rms inward flow 1-3 cells from the wall (../wallflow.npz; independent of the defects' elastic confinement) -> T ~ v_flow^2;
  (b) DEFECT: resident depth diffusivity D_r (drive_test1.csv; contaminated by confinement: a stiffer well lowers it -> cross-check only).
PASS if T_near / drive is the same across rho and coherence classes. Medians per class. -> drive_test1.txt"""
import csv, numpy as np
Z = np.load("../drag_relaxed.npz"); NS = np.array([8, 12, 18, 24, 28]); DEL = Z["DEL"]; EN = np.array([Z[f"relaxed_N{N}"] for N in NS])
def U(x, d0, rho):
    N = np.clip(rho * 200, NS[0], NS[-1]); E = np.array([np.interp(N, NS, EN[:, j]) for j in range(len(DEL))]); return np.interp(np.clip(x - d0, 0, DEL[-1]), DEL, E)
S_ = {r["run"]: r for r in csv.DictReader(open("../../shapes/shape_runs.csv"))}; DB = {r["run"]: r for r in csv.DictReader(open("drive_barrier.csv")) if not r.get("err")}
D1 = {r["run"]: r for r in csv.DictReader(open("drive_test1.csv")) if not r.get("err")}
W = np.load("../wallflow.npz"); dc = 0.5 * (W["DB"][:-1] + W["DB"][1:]); band = (dc > 1) & (dc < 3); WF = {r: a[1][band].mean() for r, a in zip(W["runs"], W["prof"])}
XG = np.array(sorted(float(k[1:]) for k in next(iter(DB.values())) if k.startswith("S"))); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); rows = []
for run, r in DB.items():
    if run not in D1 or run not in WF: continue
    s = S_[run]; rho, C, asp, n, d0 = float(s["rho"]), float(s["coh"]), float(s["aspect"]), int(r["n"]), float(r["d0"])
    S = np.array([float(r[f"S{x:g}"]) for x in XG]); ok = (XG >= d0 + 0.5) & (XG <= d0 + 2.5) & (S * n >= 5) & (S > 0)
    if ok.sum() < 3: continue
    u = np.array([U(x, d0, rho) for x in XG[ok]]); b = np.polyfit(u, np.log(S[ok]), 1)
    if b[0] >= 0: continue
    rows.append(dict(lev=int(np.argmin(abs(rho - LV))), C=C, asp=asp, Tn=-1 / b[0], vf2=WF[run] ** 2, Dr=float(D1[run]["D_r"])))
L = ["TEST 1: near-layer active temperature vs measured drive (medians per class; T_near in kBT from x in [d0+0.5, d0+2.5], relaxed-neighbour U)",
     "class | rho | runs | T_near | v_flow^2 (rms inward flow 1-3 cells)^2 | T_near / v_flow^2 (x1e-3) | D_r | T_near / D_r"]
for cname, sel in (("incoherent, aspect>=1.1", lambda q: q["C"] < 0.85 and q["asp"] > 1.05), ("coherent, aspect>=1.1", lambda q: q["C"] > 0.85 and q["asp"] > 1.05), ("circles", lambda q: q["asp"] < 1.05)):
    for i, r0 in enumerate(LV):
        g = [q for q in rows if q["lev"] == i and sel(q)]
        if len(g) < 3: continue
        md = lambda k: np.median([q[k] for q in g])
        L.append(f"{cname:24s} | {r0:.2f} | {len(g):3d} | {md('Tn'):6.1f} | {md('vf2'):.4f} | {np.median([q['Tn'] / q['vf2'] for q in g]) / 1000:6.2f} | {md('Dr'):.4f} | {np.median([q['Tn'] / q['Dr'] for q in g]):7.0f}")
open("drive_test1.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
