#!/usr/bin/env python3
"""Topological current vs width with more seeds (2026-10-16): mixed (+6/-6) and controls (+6/+6), L 192, W 20/28/40/60; 4 mixed + 2 control seeds each,
correct particle number pop = 20 L (W + 2); template channel_current/W24_mixed_s3 (6000 t.u.; mpcd_channelwind.out)."""
import json, copy, os
base = json.load(open("../channel_current/W24_mixed_s3/input.json")); L = 192; jobs = []
for W in (20, 28, 40, 60):
    for typ, ct, ns in (("mixed", 6, range(11, 15)), ("sym", -6, range(11, 13))):
        for s in ns:
            j = copy.deepcopy(base); j["domain"] = [L, W + 2]; j["species"][0]["pop"] = 20 * L * (W + 2); j["BC"][2]["Q"] = [0, 1, 0]; j["BC"][3]["Q"] = [0, 1 + W, 0]
            j["BC"][2]["customAnch"] = 6; j["BC"][3]["customAnch"] = ct; j["seed"] = 370000 + 100 * W + s + (0 if typ == "mixed" else 50)
            j["comment"] = f"CURRENT vs WIDTH (2026-10-16): L 192, W {W}, {typ}, pop {20*L*(W+2)}"
            d = f"W{W}_{typ}_s{s}"; os.makedirs(d, exist_ok=True); json.dump(j, open(d + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(d))
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs))
