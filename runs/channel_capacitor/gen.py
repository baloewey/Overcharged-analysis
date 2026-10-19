#!/usr/bin/env python3
"""Narrow mixed channels: where does the passive 'capacitor' discharge into the unscreened defect-free texture? (user, 2026-10-16)
Estimate W* ~ 14 for m 6 (unscreened ~ m^2 W, screened ~ m -> W* ~ 1/m). Mixed walls +m/-m (customAnch m on both, as channel_fix), L 192, pop 20 L (W+2),
6000 t.u. (mpcd_channelwind.out). m 6: W 11-15, passive x1 + active x2; m 3: W 20/28/36 passive; m 12: W 6/8/10 passive."""
import json, copy, os
base = json.load(open("../channel_current/W24_mixed_s3/input.json")); L = 192; jobs = []
def make(name, W, m, act, seed):
    j = copy.deepcopy(base); j["domain"] = [L, W + 2]; j["species"][0]["pop"] = 20 * L * (W + 2); j["species"][0]["act"] = act
    j["BC"][2]["Q"] = [0, 1, 0]; j["BC"][3]["Q"] = [0, 1 + W, 0]; j["BC"][2]["customAnch"] = m; j["BC"][3]["customAnch"] = m; j["seed"] = seed
    j["comment"] = f"CAPACITOR crossover (2026-10-16): mixed +{m}/-{m}, L 192, W {W}, act {act}, pop {20*L*(W+2)}"
    os.makedirs(name, exist_ok=True); json.dump(j, open(name + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(name))
for W in (11, 12, 13, 14, 15):
    make(f"m06_W{W}_passive", W, 6, 0.0, 380000 + W)
    for s in (1, 2): make(f"m06_W{W}_act_s{s}", W, 6, 0.07, 381000 + 10 * W + s)
for W in (20, 28, 36): make(f"m03_W{W}_passive", W, 3, 0.0, 382000 + W)
for W in (6, 8, 10): make(f"m12_W{W:02d}_passive", W, 12, 0.0, 383000 + W)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
