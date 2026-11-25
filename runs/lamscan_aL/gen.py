#!/usr/bin/env python3
"""Active length vs flow-alignment parameter lambda (user, 2026-11-25): periodic L200 box, zeta 0.07, shearSusc 0.5, same setup/binary as chiscan aL_chi*;
lambda (tumble) 0.8 / 1.2 / 1.5 / 2 (reference, same batch) / 3 / 4 / 6 x 2 seeds. l_a = sqrt(A / <N_defects>)."""
import json, os
base = json.load(open("../chiscan/aL_chi1.0_s1/input.json")); jobs = []
for lam in (0.8, 1.2, 1.5, 2.0, 3.0, 4.0, 6.0):
    for s in (1, 2):
        d = json.loads(json.dumps(base)); d["species"][0]["shearSusc"] = 0.5; d["species"][0]["tumble"] = lam
        d["seed"] = 79300000 + int(lam * 10) * 10 + s; d["comment"] = f"active length L200 act 0.07, shearSusc 0.5, tumble {lam} (lamscan_aL 25 Nov)"
        n = f"aL_lam{lam}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n")
