#!/usr/bin/env python3
"""Passive imposed-vortex controls for the chi scan and the U scan (2026-10-27, user). Template runs/swirl_crowd/A03_k08_s1 (circle R 31.83, swirl 0.03 within
0.6R, lambda 2, act 0, 1000 + 3000 t.u., defects every 2 t.u.). Sets: chi 0.5 / 1 at mfpot 10, and mfpot 20 at chi 0.5; k 4 / 8 / 12; 2 seeds.
Questions: does the passive crowding drop of the tilt vanish at chi 1 (as active)? Are passive near-wall excursions K-independent (mfpot 20 vs 10)?"""
import json, copy, os, random
T = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../swirl_crowd/A03_k08_s1/input.json"))); random.seed(20261027); jobs = []
for tag, chi, U in (("chi0.5_U10", 0.5, 10), ("chi1.0_U10", 1.0, 10), ("chi0.5_U20", 0.5, 20)):
    for k in (4, 8, 12):
        for s in (1, 2):
            J = copy.deepcopy(T); J["BC"][0]["customAnch"] = k; J["species"][0].update(shearSusc=chi, mfpot=U); J["seed"] = random.randrange(10**7, 10**8)
            J["comment"] = f"passive imposed vortex (swirl 0.03, 0.6R), circle R 31.83, k {k}, shearSusc {chi}, mfpot {U}; control for runs/chiscan and runs/Uscan (27 Oct 2026)"
            d = f"{tag}_k{k:02d}_s{s}"; os.makedirs(d, exist_ok=True); json.dump(J, open(d + "/input.json", "w"), indent=1); jobs.append(d)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
