#!/usr/bin/env python3
"""Stage 2, realistic active disorder (from active_wander.txt: active ellipses eps 1.5/2.0 have sigma_gap 6-8 cells and 7-8 of 10 sites occupied).
Combos: (sigma_t, vacancies) in {(6, 2), (6, 3), (8, 2), (8, 3)}, 60 samples each, eps in {1.0001, 1.2, 1.5, 2.0, 2.5}. Chunks cfg_3_cNN.json."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ns = {"__file__": os.path.join(HERE, "make_configs.py")}; exec(open(os.path.join(HERE, "make_configs.py")).read().split("out = {}")[0], ns)
wall, at = ns["wall"], ns["at"]
rng = np.random.default_rng(31337); NS, NN, DEPTH, NCHUNK = 60, 10, 1.5, 4
cfgs = []
for eps, tag in ((1.0001, "e10"), (1.2, "e12"), (1.5, "e15"), (2.0, "e20"), (2.5, "e25")):
    A = np.sqrt(2500 * eps / np.pi); B = A / eps; Wd = wall(A, B); P = Wd[2]; base = (np.arange(NN) + 0.5) * P / NN
    for st, nv in ((6.0, 2), (6.0, 3), (8.0, 2), (8.0, 3)):
        for k in range(NS):
            keep = np.sort(rng.choice(NN, NN - nv, replace=False))
            sa = np.sort((base[keep] + rng.normal(0, st, NN - nv)) % P); dp = np.clip(DEPTH + rng.normal(0, 0.5, NN - nv), 0.8, 4)
            cfgs.append(dict(id=f"real_s{st:g}_v{nv}_{k}", tag=tag, kind="realistic", sigma=st, vac=nv, A=A, B=B, sites=[at(Wd, a, d) for a, d in zip(sa, dp)]))
order = rng.permutation(len(cfgs))
for c in range(NCHUNK):
    json.dump([cfgs[i] for i in order[c::NCHUNK]], open(os.path.join(HERE, f"cfg_3_c{c:02d}.json"), "w"))
print(len(cfgs), "configs,", sum(len(c["sites"]) for c in cfgs), "solves")
