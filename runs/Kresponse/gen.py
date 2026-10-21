#!/usr/bin/env python3
"""Static response test of the Frank anisotropy (2026-10-21). Passive channel, periodic in x (L 64), no-slip walls at y = 1 and 41 (W 40),
density 20. Planar walls (director along x) + field at alpha = 10 deg from x: the director heals from 0 at the walls to ~alpha in the bulk;
gradient along y with director ~x = SPLAY, healing length xi ~ sqrt(K1 / field). Homeotropic walls (director along y) + field at 10 deg from y:
gradient along the director = BEND, xi ~ sqrt(K3 / field). Same field in both -> (xi_splay / xi_bend)^2 = K1/K3 of the RESPONSE (mean profile,
no equipartition). Coupling: production (tumble 2, shearSusc 0.5) and shearSusc 0 (control, expected isotropic). Field scan H 4, 8, 12, 16 (code divides the field by the density 20; estimate xi ~ 54/H cells)
(magnSusc 0.001 default)."""
import json, os, math
L, W, rho, alpha = 64, 40, 20, math.radians(10)
def pbc(side):
    return {"comment": f"{side} PBC", "phantom": 0, "Q": [0 if side == "Left" else L, 0, 0], "aInv": [1 if side == "Left" else -1, 0, 0],
            "P": [1, 1, 1, 1], "R": 0, "DN": L, "MVN": 1, "MVT": 1}
def wall(y, up, geom):
    mun, mut = (0, 1) if geom == "planar" else (1, 0)
    return {"comment": f"{'Bottom' if up else 'Top'} wall, {geom}", "phantom": 1, "P": [1, 1, 1, 1], "DN": 0, "MVN": -1, "MVT": -1,
            "MUN": mun, "MUT": mut, "inv": 0, "R": 0, "Q": [0, y, 0], "aInv": [0, 1 if up else -1, 0]}
jobs = []
for coup, sp in (("prod", {}), ("susc0", {"shearSusc": 0.0})):
    for geom in ("planar", "homeo"):
        for H in (4, 8, 12, 16):
            ang = alpha if geom == "planar" else math.pi / 2 - alpha
            d = {"comment": f"Frank response test: {geom} walls, field H {H} at 10 deg from the wall axis, coupling {coup} (2026-10-21)",
                 "domain": [L, W + 2], "simSteps": 15000, "warmUp": 3000, "seed": 7000 + H + (100 if geom == "homeo" else 0) + (1000 if coup == "susc0" else 0),
                 "collOp": 21, "lc": 1, "species": [dict({"pop": L * W * rho, "act": 0.0, "oDist": 1 if geom == "planar" else 2}, **sp)],
                 "BC": [pbc("Left"), pbc("Right"), wall(1, True, "planar" if geom == "planar" else "homeo"), wall(W + 1, False, "planar" if geom == "planar" else "homeo")],
                 "mag": [H * math.cos(ang), H * math.sin(ang), 0.0], "checkpointOut": 0, "synopsisOut": 1, "defectsOut": 100, "dirSOut": 20}
            n = f"{coup}_{geom}_H{H}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1)
            jobs.append(f"/home/baloewe/Overcharged/runs/Kresponse/{n}")
open("jobs.txt", "w").write("\n".join(jobs) + "\n")
