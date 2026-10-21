#!/usr/bin/env python3
"""Static response test, long-wavelength version (2026-10-21; follows runs/Kresponse). Passive channel L 64 (periodic), W 80 (walls y = 1, 81), density 20.
Planar walls + field 10 deg from x (splay healing) or homeotropic + field 10 deg from y (bend healing). Weak fields H 4, 6 (xi ~ 54/H ~ 9-14 cells,
<< half-width 39.5, so the bulk reaches the field angle and a fixed-bulk fit is clean). Production coupling (tumble 2, shearSusc 0.5), 3 seeds;
control shearSusc 0 at H 4, 2 seeds. 40000 steps after 3000 warm-up (4000 t.u.); director every 10 t.u."""
import json, os, math
L, W, rho, alpha = 64, 80, 20, math.radians(10)
def pbc(side):
    return {"comment": f"{side} PBC", "phantom": 0, "Q": [0 if side == "Left" else L, 0, 0], "aInv": [1 if side == "Left" else -1, 0, 0],
            "P": [1, 1, 1, 1], "R": 0, "DN": L, "MVN": 1, "MVT": 1}
def wall(y, up, geom):
    mun, mut = (0, 1) if geom == "planar" else (1, 0)
    return {"comment": f"{'Bottom' if up else 'Top'} wall, {geom}", "phantom": 1, "P": [1, 1, 1, 1], "DN": 0, "MVN": -1, "MVT": -1,
            "MUN": mun, "MUT": mut, "inv": 0, "R": 0, "Q": [0, y, 0], "aInv": [0, 1 if up else -1, 0]}
sets = [("prod", {}, H, s) for H in (4, 6) for s in (1, 2, 3)] + [("susc0", {"shearSusc": 0.0}, 4, s) for s in (1, 2)]
jobs = []
for coup, sp, H, s in sets:
    for geom in ("planar", "homeo"):
        ang = alpha if geom == "planar" else math.pi / 2 - alpha
        d = {"comment": f"Frank response test (long wavelength): {geom} walls, W 80, field H {H} at 10 deg from the wall axis, coupling {coup}, seed {s} (2026-10-21)",
             "domain": [L, W + 2], "simSteps": 40000, "warmUp": 3000, "seed": 9000 + 100 * s + 10 * H + (1 if geom == "homeo" else 0) + (5000 if coup == "susc0" else 0),
             "collOp": 21, "lc": 1, "species": [dict({"pop": L * W * rho, "act": 0.0, "oDist": 1 if geom == "planar" else 2}, **sp)],
             "BC": [pbc("Left"), pbc("Right"), wall(1, True, geom), wall(W + 1, False, geom)],
             "mag": [H * math.cos(ang), H * math.sin(ang), 0.0], "checkpointOut": 0, "synopsisOut": 1, "defectsOut": 1000, "dirSOut": 100}
        n = f"{coup}_{geom}_H{H}_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(f"{n}/input.json", "w"), indent=1)
        jobs.append(f"/home/baloewe/Overcharged/runs/Kresponse2/{n}")
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "runs")
