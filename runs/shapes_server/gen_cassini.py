#!/usr/bin/env python3
"""
Generate an input.json for an active nematic in a no-slip Cassini-oval cavity d1*d2 = b^2 (foci at
±a along x) with homeotropic custom anchoring of winding k. Needs the `cassini` branch build of
NIAMH-MPCD (niamh-overcharged), which adds the BC keys "shape": 1 and "cassiniA".

Same settings as runs/cavity/gen_cavity.py (no-slip, phantom, homeotropic, legacy tracker); only the
wall shape differs. Box: smallest odd integers >= 2*half-length + 2 (x) and 2*half-width + 2 (y),
oval centred on the central cell centre. Population = dens * (numerical oval area).

Geometry of the oval (b > a):
  half-length (x)  sqrt(b^2 + a^2)      waist half-width (y at x = 0)  sqrt(b^2 - a^2)
  concave neck for b/a < sqrt(2);  a = 0 is the circle of radius b.

Usage:
  gen_cassini.py OUTDIR --k 3 --b 20 --a 12.5 [--act 0] [--emode 3] [--seed 1] [--steps ...]
"""
import argparse
import json
import math
import os

import numpy as np


def shape(a, b, n=20000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    area = 0.5 * r2.mean() * 2 * np.pi
    x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    return {"area": area, "xmax": x.max(), "ymax": y.max(),
            "half_length": math.sqrt(b * b + a * a), "waist": math.sqrt(b * b - a * a)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--b", type=float, required=True, help="oval parameter (BC R)")
    ap.add_argument("--a", type=float, required=True, help="focal half-distance (cassiniA), 0 <= a < b")
    ap.add_argument("--emode", type=int, default=3, help="0 normal angle, 3 arclength by polar angle")
    ap.add_argument("--act", type=float, default=0.07)
    ap.add_argument("--dens", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--warmup", type=int, default=10000)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--defout", type=int, default=20)
    ap.add_argument("--dirout", type=int, default=500)
    ap.add_argument("--flowout", type=int, default=500)
    ap.add_argument("--checkpoint", type=int, default=10000)
    ap.add_argument("--circle", action="store_true", help="write the ordinary circle of radius b instead (reference)")
    a = ap.parse_args()
    if not 0 <= a.a < a.b:
        raise SystemExit("need 0 <= a < b")
    g = shape(a.a, a.b)
    Lx = math.ceil(2 * g["xmax"] + 2); Lx += (Lx % 2 == 0)
    Ly = math.ceil(2 * g["ymax"] + 2); Ly += (Ly % 2 == 0)
    if a.circle:
        Lx = Ly = max(Lx, Ly)
    bc = {
        "comment":    "No-slip Cassini cavity, homeotropic custom anchoring" if not a.circle else "No-slip circular cavity (reference)",
        "phantom":    1,
        "Q":          [Lx / 2, Ly / 2, 0],
        "aInv":       [1, 1, 1],
        "P":          [2, 2, 2, 2],
        "R":          round(a.b, 6),
        "DN":         0,
        "MVN":        -1,
        "MVT":        -1,
        "MUN":        1,
        "MUT":        0,
        "inv":        1,
        "customAnch": a.k,
    }
    if not a.circle:
        bc.update({"shape": 1, "cassiniA": round(a.a, 6), "eMode": a.emode})
    area = g["area"] if not a.circle else math.pi * a.b ** 2
    inp = {
        "comment":       (f"cassini k={a.k} b={a.b} a={a.a} b/a={a.b / a.a if a.a else float('inf'):.3f} "
                          f"eMode={a.emode} act={a.act}") if not a.circle else f"circle k={a.k} R={a.b} act={a.act}",
        "domain":        [Lx, Ly],
        "simSteps":      a.steps,
        "warmUp":        a.warmup,
        "seed":          a.seed,
        "collOp":        21,
        "lc":            1,
        "species":       [{"pop": round(a.dens * area), "act": a.act}],
        "BC":            [bc],
        "checkpointOut": a.checkpoint,
        "synopsisOut":   1,
        "defectTracker": 0,
    }
    for key, val in (("defectsOut", a.defout), ("dirSOut", a.dirout), ("flowOut", a.flowout)):
        if val > 0:
            inp[key] = val
    os.makedirs(a.outdir, exist_ok=True)
    json.dump(inp, open(os.path.join(a.outdir, "input.json"), "w"), indent=4)
    print(f"{a.outdir}: domain {Lx}x{Ly}, area {area:.1f}, pop {inp['species'][0]['pop']}, "
          f"half-length {g['half_length']:.2f}, waist {g['waist']:.2f}")


if __name__ == "__main__":
    main()
