#!/usr/bin/env python3
"""Layer stationarity by wall density and aspect (shape_fourier.csv): f_lat = static fraction at the lattice mode 2k-2 (1 = pinned lattice), S_lat, decay time.
Medians over runs (R_eq and seeds pooled). -> analyse_fourier.txt"""
import csv, numpy as np
F = {r["run"]: r for r in csv.DictReader(open("shape_fourier.csv")) if not r.get("err")}; S = {r["run"]: r for r in csv.DictReader(open("../../shapes/shape_runs.csv"))}
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); CL = (("circle", 1.0, 1.0), ("aspect 1.1-1.3", 1.05, 1.35), ("aspect 1.35-1.7", 1.35, 1.75), ("aspect 1.9-2.5", 1.85, 3))
L = ["LAYER STATIONARITY (wall Fourier at the lattice mode m = 2k-2): f_lat median [S_lat median] (n runs); columns rho"]
L.append("class            | " + " | ".join(f"rho {r:.2f}        " for r in LV))
for nm, lo, hi in CL:
    row = []
    for i, r0 in enumerate(LV):
        rs = [r for r in F if lo <= float(S[r]["aspect"]) <= hi and np.argmin(abs(float(S[r]["rho"]) - LV)) == i]
        f = np.array([float(F[r]["f_lat"]) for r in rs]); s = np.array([float(F[r]["S_lat"]) for r in rs])
        row.append(f"{np.median(f):.2f} [{np.median(s):4.1f}] ({len(rs):2d})" if len(rs) else "-")
    L.append(f"{nm:16s} | " + " | ".join(row))
open("analyse_fourier.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
