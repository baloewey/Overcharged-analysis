#!/usr/bin/env python3
"""Pressure-driven winding channels (user, 2026-10-16: impose a pressure difference along the channel). A uniform body force g along x ('grav') is the
periodic-channel equivalent of a pressure gradient dp/dx = -rho g. Questions: does the flow respond differently to +g and -g (rectification: the mixed
channel's walls wind the director the same way in the lab frame, so x -> -x is not a symmetry), does the forcing change the wall layers, the release of
charge and the topological current? L 192, W 24, pop 20 L (W+2), 6000 t.u. (mpcd_channelwind.out). Sets: mixed (+6/-6) and symmetric (+6/+6) active at
g = +-0.001 and +-0.003 (2 seeds each); passive mixed and homeotropic (m 0) active at g = +-0.003 (1 seed each) as controls."""
import json, copy, os
base = json.load(open("../channel_current/W24_mixed_s3/input.json")); L, W = 192, 24; jobs = []
def make(name, cb, ct_anch, act, g, seed):
    j = copy.deepcopy(base); j["domain"] = [L, W + 2]; j["species"][0]["pop"] = 20 * L * (W + 2); j["species"][0]["act"] = act
    j["BC"][2]["customAnch"] = cb; j["BC"][3]["customAnch"] = ct_anch; j["grav"] = [g, 0.0, 0.0]; j["seed"] = seed
    j["comment"] = f"PRESSURE-DRIVEN channel (2026-10-16): L 192, W 24, walls {cb:+d}/{-ct_anch:+d}, act {act}, body force g_x {g:+.4f}"
    os.makedirs(name, exist_ok=True); json.dump(j, open(name + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(name))
i = 0
for typ, ct in (("mixed", 6), ("sym", -6)):
    for g in (0.001, -0.001, 0.003, -0.003):
        for s in (1, 2): i += 1; make(f"{typ}_g{g:+.3f}_s{s}", 6, ct, 0.07, g, 390000 + i)
for g in (0.003, -0.003):
    i += 1; make(f"mixed_passive_g{g:+.3f}", 6, 6, 0.0, g, 390000 + i)
    i += 1; make(f"homeo_g{g:+.3f}", 0, 0, 0.07, g, 390000 + i)
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
