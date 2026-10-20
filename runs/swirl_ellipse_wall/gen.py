#!/usr/bin/env python3
"""Wall-reaching imposed swirl in ACTIVE ellipses (user, 2026-10-20). The earlier swirl_ellipse(_strong) runs applied the azimuthal body force only inside a disk of
radius 0.6 x the minor semi-axis (12-14 cells), so the wall flow stayed at 0.09-0.17 of v_rms and the layer was never repaired. Here swirlRadius = 0 (force everywhere,
as the uniform body force that re-screens active channels). Target: near-wall flow ~0.5 v_rms as in active circles (0.38/0.74) and forced channels.
Templates: swirl_ellipse_strong/active_eps{1.5,2.0}_A0.06 (k 6, area 2500, act 0.07, lambda 2, 4000 t.u. after 1000 t.u. warm-up, flow every 20 t.u.). Binary bin/mpcd_swirl.out."""
import json, copy, os
jobs = []
for eps in ("1.5", "2.0"):
    base = json.load(open(f"../swirl_ellipse_strong/active_eps{eps}_A0.06/input.json"))
    for A in (0.0005, 0.001):
        for s in (1, 2):
            j = copy.deepcopy(base); j["swirlAcc"] = A; j["swirlRadius"] = 0.0; j["seed"] = 820000 + int(float(eps) * 100) * 100 + int(A * 10000) * 10 + s
            j["comment"] = f"WALL-REACHING SWIRL (2026-10-20): active ellipse eps {eps}, k 6, area 2500, swirlAcc {A} everywhere (swirlRadius 0), seed s{s}"
            d = f"active_eps{eps}_A{A}_s{s}"; os.makedirs(d, exist_ok=True); json.dump(j, open(d + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(d))
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
