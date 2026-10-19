#!/usr/bin/env python3
"""Channel re-runs with the CORRECT particle number (2026-10-16). channel_current and channel_crowd copied pop 99840 (right only for W 24) to every width:
fluid density 32.5 (W 16), 16.2 (W 32), 13.0 (W 40) instead of ~21. Standard: pop = 20 L (W + 2). Also the crowded channels with the walls moved 4 cells away
from the box edges (charge did not close for m >= 18 at the bottom wall). Template: channel_current/W24_mixed_s3 (mpcd_channelwind.out)."""
import json, copy, os
base = json.load(open("../channel_current/W24_mixed_s3/input.json")); L = 192; jobs = []
def make(name, W, cb, ct_anch, seed, off=1, comment=""):
    j = copy.deepcopy(base); j["domain"] = [L, W + 2 * off]; j["species"][0]["pop"] = 20 * L * (W + 2)
    j["BC"][2]["Q"] = [0, off, 0]; j["BC"][3]["Q"] = [0, off + W, 0]; j["BC"][2]["customAnch"] = cb; j["BC"][3]["customAnch"] = ct_anch; j["seed"] = seed
    j["comment"] = comment + f" (pop 20 L (W+2) = {20*L*(W+2)}; walls at y {off} and {off+W})"
    os.makedirs(name, exist_ok=True); json.dump(j, open(name + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(name))
# (a) topological current at W 16 and W 32 with the correct density
for W, ns_m, ns_s in ((16, range(3, 9), range(3, 7)), (32, range(1, 5), range(1, 3))):
    for s in ns_m: make(f"W{W}_mixed_s{s}", W, 6, 6, 330000 + 100 * W + s, comment=f"TOPOLOGICAL CURRENT re-run: L 192, W {W}, wall charges +6/-6")
    for s in ns_s: make(f"W{W}_sym_s{s}", W, 6, -6, 340000 + 100 * W + s, comment=f"control re-run: L 192, W {W}, wall charges +6/+6")
# (b) crowding, walls 4 cells from the box edges; W 40 also corrected density
for W, ms in ((24, (18, 24)), (40, (12, 18, 24))):
    for m in ms:
        for s in (1, 2): make(f"crowd_W{W}_m{m}_mixed_s{s}", W, m, m, 350000 + 1000 * m + 10 * W + s, off=4, comment=f"CROWDING re-run: L 192, W {W}, +{m}/-{m}")
for m in (18, 24):
    for s in (1, 2): make(f"crowd_W40_m{m}_sym_s{s}", 40, m, -m, 360000 + 1000 * m + s, off=4, comment=f"CROWDING re-run: L 192, W 40, +{m}/+{m}")
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
