#!/usr/bin/env python3
"""Intro-slide runs for the workshop talk (2026-11-23, user: softer intro for non-LC active-matter audience). Director + defects every 2 t.u.
  bulk80     periodic 80 x 80 box, act 0.07: active turbulence with +-1/2 defects (left panel)
  k1_active  circle R 33.315 (R_eq 1.5 l_a), customAnch 1 (homeotropic, winding 1), act 0.07, 2000 t.u.: prediction 2k-2 = 0 wall defects, bulk +1
  k1_passive same circle, act 0: the passive reference
Templates: shapes/circle/circle/Req1.5/k06/s01/input.json (circle), Uscan/aL_U10_s1/input.json (bulk)."""
import json, os
C = json.load(open(os.path.expanduser("~/Overcharged/runs/shapes/circle/circle/Req1.5/k06/s01/input.json")))
B = json.load(open(os.path.expanduser("~/Overcharged/runs/Uscan/aL_U10_s1/input.json")))
jobs = {}
jobs["bulk80"] = dict(B, domain=[80, 80], warmUp=10000, simSteps=3000, seed=83000001, dirSOut=20, flowOut=20, defectsOut=20, checkpointOut=0,
                      comment="talk intro: periodic 80x80 active turbulence, act 0.07, director every 2 t.u.")
for name, act, steps, seed, dout in (("k1_active", 0.07, 20000, 83000002, 20), ("k1_passive", 0.0, 3000, 83000003, 200)):
    d = json.loads(json.dumps(C)); d["species"][0]["act"] = act; d["BC"][0]["customAnch"] = 1; d["BC"][0]["comment"] = "No-slip circle, homeotropic winding 1"
    d.update(seed=seed, simSteps=steps, dirSOut=dout, flowOut=20, defectsOut=20, checkpointOut=0, comment=f"talk intro: circle R 33.315 k=1 act={act}")
    jobs[name] = d
for n, d in jobs.items():
    os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=4)
open("jobs.txt", "w").write("\n".join(jobs) + "\n")
