#!/usr/bin/env python3
"""Test 2 (user, 2026-11-25): passive MPCD layers under an imposed vortex whose SHEAR is confined to a thin layer at the wall (annular forcing,
swirlInnerRadius; bin/mpcd_swirlring.out, branch swirlring). Inside the annulus the fluid rotates rigidly (no strain: no flow-aligned band behind the layer);
all the shear sits between the annulus and the no-slip wall, where the flow is simple shear (aligned angle theta_L for any chi).
Total torque matched to the disk forcing of runs/swirl_crowd A03 (swirl 0.03 within 0.6R): A' = 0.03 rs^3 / (r2^3 - r1^3), so the wall stress is the same.
  thin:   annulus R-5 .. R-3 (shear in the outer 3 cells), A' = 0.045
  medium: annulus R-9 .. R-7 (shear in the outer 7 cells), A' = 0.061
  disk:   control as A03 (0.03 within 0.6R) for k 4 (A03 exists for k 6-16)
k 4/6/8/10/13/16, 2 seeds; R 31.83, lambda 2, act 0, 1000 + 3000 t.u.; defects every 2 t.u., director every 20 t.u., flow every 10 t.u.
Predictions. Picture (sparse set by the near-wall simple shear, crowded by the band): sparse tilts ~2 theta_L in all conditions. In thin/medium there is no
strained band behind the layer, so crowded layers, whose near field (~2 cells) sits inside the shear layer, should also stay near 2 theta_L (drop reduced or
gone). Alternative (sparse overtilt = the defect's own rotated neighbourhood; drop intrinsic to crowding): sparse overtilt and drop persist."""
import json, os
T = json.load(open("/home/baloewe/Overcharged/runs/swirl_crowd/A03_template.json")); R = T["BC"][0]["R"]; rs = T["swirlRadius"]
jobs = []
def make(name, k, s, acc, rin, rout, note):
    d = json.loads(json.dumps(T)); d["BC"][0]["customAnch"] = k
    d.update(seed=95000000 + k * 100 + s + {"thin": 10, "medium": 20, "disk": 30}[name] * 1000, swirlAcc=acc, swirlRadius=rout, swirlInnerRadius=rin, dirSOut=200,
             comment=f"swirl_thin (test 2, 2026-11-25): PASSIVE circle R {R:.2f}, k {k}, {note}")
    n = f"{name}_k{k:02d}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(n + "/input.json", "w"), indent=4); jobs.append(n)
for k in (4, 6, 8, 10, 13, 16):
    for s in (1, 2):
        for name, r1, r2 in (("thin", R - 5, R - 3), ("medium", R - 9, R - 7)):
            acc = round(0.03 * rs ** 3 / (r2 ** 3 - r1 ** 3), 4); make(name, k, s, acc, r1, r2, f"annular swirl {acc} for {r1:.2f} <= r <= {r2:.2f} (shear confined to the outer {R - r2:.0f} cells)")
        if k == 4: make("disk", k, s, 0.03, 0.0, rs, "disk swirl 0.03 within 0.6R (A03 control)")
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
