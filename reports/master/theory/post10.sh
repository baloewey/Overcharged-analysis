#!/bin/bash
# Second stage of the 10-seed reanalysis (27 Oct): merge server flow results, summaries, figures. Run after rerun10.sh and the server flow job.
set -e
T=/home/baloewe/Physics/Overcharged/reports/master/theory; cd $T
scp -q baloewe@146.155.46.150:~/Overcharged/analysis_flow10/flow10.csv shapes/flow10/ && scp -q baloewe@146.155.46.150:~/Overcharged/analysis_flow10/wallflow10.npz shapes/flow10/
python3 - <<'PY'
import csv, numpy as np
F = {r["run"]: r for r in csv.DictReader(open("shapes/flow10/flow10.csv"))}; R = list(csv.DictReader(open("shapes/shape_runs.csv"))); cols = list(R[0].keys())
for r in R:
    f = F.get(r["run"])
    if f and not f["err"]: r["circ"], r["two"] = f["circ"], f["two"]
with open("shapes/shape_runs.csv", "w", newline="") as fo: w = csv.DictWriter(fo, fieldnames=cols); w.writeheader(); w.writerows(R)
print("merged flow for", sum(1 for r in R if r["circ"] not in ("", "nan")), "runs")
PY
cp --update=none crowd_protect/wallflow.npz crowd_protect/wallflow_3seeds.npz; cp shapes/flow10/wallflow10.npz crowd_protect/wallflow.npz
