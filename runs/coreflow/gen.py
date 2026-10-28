#!/usr/bin/env python3
"""Flow-at-core test, MPCD (2026-10-28, user). Same runs as the shape sweep aspect-2 ellipses R_eq 1.5 l_a, k 5 (sparse) and 14 (crowded), act 0.07
(templates runs/shapes/ellipse/eps2.00/Req1.5/k05|k14/s01/input.json), but flowOut 20 (flow averaged over 2 t.u., the defect cadence) and simSteps 50000
(5000 t.u. after the 10000-step warm-up), dirSOut 1000. 5 new seeds per k. Question: at wall cores, does the inward-normal flow (the drive) or the
defect velocity relative to the flow (elastic + self-propulsion) change with crowding?"""
import json, os
T = os.path.expanduser("~/Overcharged/runs/shapes/ellipse/eps2.00/Req1.5")
jobs = []
for k in (5, 14):
    d = json.load(open(f"{T}/k{k:02d}/s01/input.json"))
    for s in range(1, 6):
        n = f"ell2_k{k:02d}_s{s}"; os.makedirs(n, exist_ok=True)
        e = dict(d, seed=82000000 + 100 * k + s, flowOut=20, simSteps=50000, comment=d["comment"] + " | coreflow flowOut 20")
        json.dump(e, open(f"{n}/input.json", "w"), indent=4); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n")
