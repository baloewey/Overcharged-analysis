#!/usr/bin/env python3
"""Passive -> active switch (user, 2026-10-21): watch the tilt-position lock and the registry form when activity is turned on.
Start from the final checkpoints (t = 3000 t.u.) of the 8 passive k = 6, R = 31.83 runs of runs/pinning_ic (different initial directors, so
different frozen passive registries, 12 wall defects). In each checkpoint: activity 0 -> 0.07, total steps 30000 -> 70000 (4000 t.u. active),
flow output every 100 steps (10 t.u.); defects every 2 t.u., director every 100 t.u. as before. Prediction (lock + flow selection, no anisotropy
needed): two defects leave, the 10-ring tilts to ~+-63 deg and every run ends at the same registry (arg Z_10 = +-1.17 rad, D = 0), as in
runs/pinning_ic_active, which was active from the start."""
import json, os
SRC = "/home/baloewe/Physics/Overcharged/runs/pinning_ic"; jobs = []
for ic in sorted(os.listdir(SRC)):
    L = open(f"{SRC}/{ic}/checkpoint.dat").read().split("\n"); assert L[0].strip() == "30000"; L[0] = "70000"
    i = [k for k, l in enumerate(L) if l.startswith("1.000000 63662 ")]; assert len(i) == 1; f = L[i[0]].split(); assert f[10] == "0.000000"; f[10] = "0.070000"; L[i[0]] = " ".join(f)
    j = [k for k, l in enumerate(L) if l.startswith("3 0 0 0 0 0 0 0 0 0 1000")]; assert len(j) == 1; g = L[j[0]].split(); g[4] = "100"; L[j[0]] = " ".join(g)
    os.makedirs(ic, exist_ok=True); open(f"{ic}/checkpoint_in.dat", "w").write("\n".join(L))
    d = json.load(open(f"{SRC}/{ic}/input.json")); d.update({"checkpointIn": "checkpoint_in.dat", "seed": -1, "checkpointOut": 0,
        "comment": f"passive -> active switch at t = 3000 from runs/pinning_ic/{ic} (act 0.07 set in the checkpoint), 2026-10-21"})
    json.dump(d, open(f"{ic}/input.json", "w"), indent=1); jobs.append(f"/home/baloewe/Overcharged/runs/switchon/{ic}")
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
