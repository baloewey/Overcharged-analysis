#!/usr/bin/env python3
"""chi = 1 circles at matched R/l_a (user, 2026-11-23). l_a(chi 1) = 15.1 (aL_chi1.0_s1/s2: 15.12, 15.09; merged+neutral estimator) vs 22.2 at chi 0.5.
Same R/l_a and k as the chi 0.5 main sweep -> same wall spacing in units of l_a; resolution limit rho = k/(pi R) <= 0.15.
  R/l_a 1.43 (R 21.6): k 4/6/8/10 x 4 seeds | R/l_a 2.0 (R 30.2): k 4/6/8/10/12 x 3 | R/l_a 1.0 (R 15.1): k 3/4/6 x 3
Template: runs/chiscan/circ_k04_chi1.0_s1/input.json (shearSusc 1.0, act 0.07, 100000 + 10000 steps, defects every 2 t.u., director + flow every 20 t.u.).
domain = ceil(2R) + 3, centre at domain/2, pop = 20 pi R^2. Question: do the overtilt, the crowding drop and the size dependence of the tilt survive at chi 1?"""
import json, math, os
T = json.load(open(os.path.expanduser("~/Overcharged/runs/chiscan/circ_k04_chi1.0_s1/input.json"))); LA = 15.1; jobs = []
for rl, ks, ns in ((1.43, (4, 6, 8, 10), 4), (2.0, (4, 6, 8, 10, 12), 3), (1.0, (3, 4, 6), 3)):
    R = round(rl * LA, 4); dom = math.ceil(2 * R) + 3
    for k in ks:
        assert k / (math.pi * R) <= 0.15 + 1e-9, (R, k)
        for s in range(1, ns + 1):
            d = json.loads(json.dumps(T)); n = f"R{rl:.2f}_k{k:02d}_s{s}"
            d.update(domain=[dom, dom], seed=84000000 + int(rl * 100) * 1000 + k * 10 + s, comment=f"chi1_matched: circle R={R} (R/l_a {rl}, l_a 15.1) k={k} act=0.07 shearSusc 1.0")
            d["species"][0]["pop"] = int(round(20 * math.pi * R * R)); d["BC"][0].update(R=R, Q=[dom / 2, dom / 2, 0], customAnch=k)
            os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=4); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
