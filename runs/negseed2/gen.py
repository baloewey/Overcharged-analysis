#!/usr/bin/env python3
"""Pair-creation-cost test (user, 2026-10-19). Seeded starts (branch seeddir, oDist 10) in mirror-sweep cavities where the -1/2 layer does / does not gain:
  sparse_circle   circle R_eq 2.0 k -6  (rho 0.043; sweep gain 0.27)
  sparse_ellipse  ellipse eps 2 R_eq 2.0 k -6 (sweep gain -0.41)
  small_circle    circle R_eq 1.0 k -10 (rho 0.143; sweep gain 1.32, saturated)
  crowded_circle  circle R_eq 2.0 k -13 (rho 0.093; sweep gain 2.03, control)
Starts: 'base' 2|k| wall -1/2, bulk 0; 'full' 2|k|+2 wall -1/2 + 2 bulk +1/2 (bulk +1). 2 seeds, act 0.07, 6000 t.u.
Barrier (pair-creation cost) hypothesis: 'full' stays full in the sparse/small cavities while 'base' stays low (history dependence).
Steady-state (not a barrier): both starts converge to the sweep value."""
import json, copy, os, subprocess
SRC = "/home/baloewe/Physics/Overcharged/runs/shapes_mirror/"
CASES = {"sparse_circle": "circle/circle/Req2.0/km06/s01", "sparse_ellipse": "ellipse/eps2.00/Req2.0/km06/s01",
         "small_circle": "circle/circle/Req1.0/km10/s01", "crowded_circle": "circle/circle/Req2.0/km13/s01"}
jobs = []
for name, src in CASES.items():
    base = json.load(open(SRC + src + "/input.json")); k = base["BC"][0]["customAnch"]
    for start in ("base", "full"):
        for s in (1, 2):
            j = copy.deepcopy(base); j["simSteps"] = 60000; j["warmUp"] = 0; j["checkpointOut"] = 0; j["seed"] = 96000000 + 1000 * abs(k) + 100 * (start == "full") + 10 * len(jobs) + s
            j["comment"] = f"PAIR-COST TEST (2026-10-19): {name} from mirror sweep {src}, start {start}, seeded, act 0.07, 6000 t.u."
            d = f"{name}_{start}_s{s}"; os.makedirs(d, exist_ok=True); json.dump(j, open(d + "/input.json", "w"), indent=1)
            args = ["--wall", str(2 * abs(k) + (2 if start == "full" else 0))] + (["--nplus", "2"] if start == "full" else [])
            subprocess.run(["python3", "/home/baloewe/Physics/Overcharged/reports/master/theory/seeding/seed_director.py", d] + args, check=True, capture_output=True)
            jobs.append(os.path.abspath(d))
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
