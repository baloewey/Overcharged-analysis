#!/usr/bin/env python3
"""chi-scan (27 Oct; tests whether the crowding drop of the tilt follows d/l_E, l_E ~ 1/sqrt(chi_HI)): shearSusc 0.25, 0.5 (production), 1.
(a) circles R 31.83 (1.43 l_a; as cavity/server_main template) at k 4, 8, 12 (rho 0.04, 0.08, 0.12), 3 seeds; director and flow every 20 t.u.
(b) active length L200 act 0.07 at shearSusc 0.25 and 1, 2 seeds."""
import json, os, copy
R = "/home/baloewe/Physics/Overcharged/runs"; jobs = []
T = json.load(open(f"{R}/cavity/server_main/main_act0.07/rho0.06/k06/s01/input.json"))
for chi in (0.25, 0.5, 1.0):
    for k in (4, 8, 12):
        for s in (1, 2, 3):
            d = copy.deepcopy(T); d["BC"][0]["customAnch"] = k; d["species"][0]["shearSusc"] = chi; d["seed"] = 78000000 + int(chi * 100) * 1000 + k * 10 + s
            d["dirSOut"] = 200; d["flowOut"] = 200; d.pop("enstrophyOut", None); d["comment"] = f"cavity k={k} R=31.8310 act=0.07, chiscan shearSusc {chi} (27 Oct 2026)"
            n = f"circ_k{k:02d}_chi{chi}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
A = json.load(open(f"{R}/activeLength/L200/act_0.070/input.json"))
for chi in (0.25, 1.0):
    for s in (1, 2):
        d = copy.deepcopy(A); d["species"][0]["shearSusc"] = chi; d["seed"] = 78200000 + int(chi * 100) * 10 + s; d["checkpointOut"] = 0; d["comment"] = f"active length L200 act 0.07, shearSusc {chi} (chiscan 27 Oct)"
        n = f"aL_chi{chi}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
