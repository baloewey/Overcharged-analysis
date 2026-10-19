#!/usr/bin/env python3
"""Activity robustness of the main claims (2026-10-16, user): circles k 6 at R = 1.43 l_a and aspect-2 ellipses at R_eq = 1.5 l_a, k 7 (wall-law winding
q ~ 0.03), at zeta 0.06 and 0.08 (l_a = 22.21 sqrt(0.07/zeta); sizes scaled so R/l_a is as at 0.07), 3 seeds, 6000 t.u. after 1000 warm-up.
Templates: server_main main_act0.07/rho0.06/k06 (circle) and shapes_server ellipse eps2.00/Req1.5/k07 (pop = 20 x area). Binary mpcd_wallact2.out."""
import json, copy, os, math
jobs = []
for zeta in (0.06, 0.08):
    f = math.sqrt(0.07 / zeta)
    for shape, tmpl in (("circle", "../cavity/server_main/main_act0.07/rho0.06/k06/s01/input.json"), ("ellipse", "../shapes_server/ellipse/eps2.00/Req1.5/k07/s01/input.json")):
        base = json.load(open(tmpl))
        for s in (1, 2, 3):
            j = copy.deepcopy(base); bc = j["BC"][0]; R = bc["R"] * f; ax = bc["aInv"][0]
            dx, dy = int(math.ceil(2 * R * ax + 3)), int(math.ceil(2 * R + 3)); j["domain"] = [dx, dy]; bc["R"] = round(R, 6); bc["Q"] = [dx / 2, dy / 2, 0]
            j["species"][0]["act"] = zeta; j["species"][0]["pop"] = int(round(20 * math.pi * R * R * ax)); j["simSteps"] = 60000; j["warmUp"] = 10000
            j.pop("enstrophyOut", None); j["seed"] = 76000000 + int(zeta * 1000) * 100 + (0 if shape == "circle" else 50) + s
            j["comment"] = f"ACTIVITY ROBUSTNESS (2026-10-16): {shape}, zeta {zeta}, sizes x{f:.4f} so R/l_a as at 0.07 (l_a {22.21*f:.2f}); k {bc['customAnch']}"
            d = f"{shape}_z{zeta:.2f}_s{s}"; os.makedirs(d, exist_ok=True); json.dump(j, open(d + "/input.json", "w"), indent=1); jobs.append(os.path.abspath(d))
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs))
