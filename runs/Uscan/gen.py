#!/usr/bin/env python3
"""U-scan (27 Oct; tests the elastic barrier of crowding protection): nematic mean-field potential mfpot 5 and 20 (production 10).
(a) cavities: aspect-2 ellipses, R_eq 1.5 l_a (as shapes_server/ellipse/eps2.00/Req1.5), k 5 (rho 0.04) and k 14 (rho 0.12), 3 seeds; the mfpot 10 references
    are the 10-seed shape sweep runs. (b) passive Frank calibration L64 (as runs/Kcalib/L64), 2 seeds. (c) active length L200 act 0.07 (as runs/activeLength/L200), 2 seeds."""
import json, os, copy
R = "/home/baloewe/Physics/Overcharged/runs"; jobs = []
for U in (5, 20):
    for k in (5, 14):
        T = json.load(open(f"{R}/shapes_server/ellipse/eps2.00/Req1.5/k{k:02d}/s01/input.json"))
        for s in (1, 2, 3):
            d = copy.deepcopy(T); d["species"][0]["mfpot"] = U; d["seed"] = 77000000 + U * 1000 + k * 10 + s
            d["comment"] = T["comment"] + f"; Uscan mfpot {U} (27 Oct 2026)"; n = f"ell2_k{k:02d}_U{U}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
    T = json.load(open(f"{R}/Kcalib/L64_s1/input.json"))
    for s in (1, 2):
        d = copy.deepcopy(T); d["species"][0]["mfpot"] = U; d["seed"] = 77100000 + U * 10 + s; d["comment"] = f"passive Frank calibration L64, mfpot {U} (Uscan 27 Oct)"
        n = f"Kcalib_U{U}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
    T = json.load(open(f"{R}/activeLength/L200/act_0.070/input.json"))
    for s in (1, 2):
        d = copy.deepcopy(T); d["species"][0]["mfpot"] = U; d["seed"] = 77200000 + U * 10 + s; d["checkpointOut"] = 0; d["comment"] = f"active length L200 act 0.07, mfpot {U} (Uscan 27 Oct)"
        n = f"aL_U{U}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
