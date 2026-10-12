#!/usr/bin/env python3
"""
Generate an input.json for an active nematic in a no-slip circular cavity with
homeotropic custom anchoring of winding k (customAnch hack).

Geometry: square box of odd side L = smallest odd integer >= 2R+2, cavity
centred on the central cell centre (L/2, L/2). Population = dens * pi R^2
(set explicitly: the "dens" key would fill the whole square box).

Either give R directly, or give rho (bare linear defect density = 2|k|/(2 pi R))
and R is set to |k|/(pi rho). Negative k (mirror runs) forces -1/2 defects.

Usage:
    gen_cavity.py OUTDIR --k 5 --rho 0.1 [--act 0.07] [--seed 1] [--dirout 500]
    gen_cavity.py OUTDIR --k 5 --R 31.83 ...
"""
import argparse
import json
import math
import os

R_MIN = 5.0


def geometry(k, rho=None, R=None, dens=20.0):
    """Return R, rho, L, centre C and population for a cavity."""
    # rho uses |k|: negative winding forces 2|k| defects just like positive
    if R is None:
        R = abs(k) / (math.pi * rho)
    rho = abs(k) / (math.pi * R)
    L = math.ceil(2 * R + 2)
    L += (L % 2 == 0)
    return {"R": R, "rho": rho, "L": L, "C": L / 2, "pop": round(dens * math.pi * R * R)}


def make_input(k, rho=None, R=None, act=0.07, dens=20.0, seed=1, warmup=10000,
               steps=100000, defout=50, dirout=0, flowout=0, enstout=0,
               checkpoint=10000):
    g = geometry(k, rho, R, dens)
    if g["R"] < R_MIN:
        raise ValueError(f"R = {g['R']:.2f} < {R_MIN}: too small to resolve")
    C = g["C"]
    inp = {
        "comment":      f"cavity k={k} R={g['R']:.4f} rho={g['rho']:.4f} act={act}",
        "domain":       [g["L"], g["L"]],
        "simSteps":     steps,
        "warmUp":       warmup,
        "seed":         seed,
        "collOp":       21,
        "lc":           1,
        "species":      [{"pop": g["pop"], "act": act}],
        "BC": [{
            "comment":    "No-slip circular cavity, homeotropic custom anchoring",
            "phantom":    1,
            "Q":          [C, C, 0],
            "aInv":       [1, 1, 1],
            "P":          [2, 2, 2, 2],
            "R":          round(g["R"], 6),
            "DN":         0,
            "MVN":        -1,
            "MVT":        -1,
            "MUN":        1,
            "MUT":        0,
            "inv":        1,
            "customAnch": k,
        }],
        "checkpointOut": checkpoint,
        "synopsisOut":   1,
    }
    for key, val in (("defectsOut", defout), ("dirSOut", dirout),
                     ("flowOut", flowout), ("enstrophyOut", enstout)):
        if val > 0:
            inp[key] = val
    return inp, g


def write_input(outdir, inp):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "input.json"), "w") as f:
        json.dump(inp, f, indent=4)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--k", type=int, required=True, help="anchoring winding (customAnch); negative for mirror runs")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--rho", type=float, help="bare linear defect density 2k/(2 pi R)")
    g.add_argument("--R", type=float, help="cavity radius")
    ap.add_argument("--act", type=float, default=0.07)
    ap.add_argument("--dens", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--warmup", type=int, default=10000)
    ap.add_argument("--steps", type=int, default=100000)
    ap.add_argument("--defout", type=int, default=50, help="defectsOut interval (0=off)")
    ap.add_argument("--dirout", type=int, default=0, help="dirSOut interval (0=off)")
    ap.add_argument("--flowout", type=int, default=0, help="flowOut interval (0=off)")
    ap.add_argument("--enstout", type=int, default=0, help="enstrophyOut interval (0=off)")
    ap.add_argument("--checkpoint", type=int, default=10000)
    a = ap.parse_args()

    try:
        inp, geo = make_input(a.k, a.rho, a.R, a.act, a.dens, a.seed, a.warmup, a.steps,
                              a.defout, a.dirout, a.flowout, a.enstout, a.checkpoint)
    except ValueError as e:
        raise SystemExit(str(e))
    write_input(a.outdir, inp)
    print(f"{a.outdir}: k={a.k} R={geo['R']:.3f} rho={geo['rho']:.4f} "
          f"L={geo['L']} pop={geo['pop']}")


if __name__ == "__main__":
    main()
