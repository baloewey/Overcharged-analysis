#!/usr/bin/env python3
"""Ellipse cavity input (semi-axes A = eps*B along x, B along y; aInv = [eps, 1, 1], R = B), eMode 2 (arclength),
homeotropic custom anchoring k; other settings as gen_cassini.py. Usage: gen_ellipse.py OUTDIR --area 2500 --eps 1.5 --k 6 ..."""
import argparse, json, math, os
ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("--area", type=float, required=True)
ap.add_argument("--eps", type=float, required=True); ap.add_argument("--k", type=int, required=True)
ap.add_argument("--act", type=float, default=0.07); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--warmup", type=int, default=10000); ap.add_argument("--steps", type=int, default=100000)
ap.add_argument("--defout", type=int, default=20); ap.add_argument("--dirout", type=int, default=1000); ap.add_argument("--flowout", type=int, default=1000)
a = ap.parse_args()
A = math.sqrt(a.area * a.eps / math.pi); B = A / a.eps
if B * B / A < 5: raise SystemExit(f"tip radius {B*B/A:.1f} < 5")
Lx = math.ceil(2 * A + 2); Lx += (Lx % 2 == 0); Ly = math.ceil(2 * B + 2); Ly += (Ly % 2 == 0)
inp = {"comment": f"ellipse k={a.k} A={A:.4f} B={B:.4f} eps={a.eps} eMode=2 act={a.act}", "domain": [Lx, Ly],
       "simSteps": a.steps, "warmUp": a.warmup, "seed": a.seed, "collOp": 21, "lc": 1,
       "species": [{"pop": round(20 * math.pi * A * B), "act": a.act}],
       "BC": [{"comment": "No-slip elliptical cavity, homeotropic custom anchoring", "phantom": 1, "Q": [Lx / 2, Ly / 2, 0],
               "aInv": [a.eps, 1, 1], "P": [2, 2, 2, 2], "R": round(B, 6), "DN": 0, "MVN": -1, "MVT": -1, "MUN": 1, "MUT": 0,
               "inv": 1, "customAnch": a.k, "eMode": 2}],
       "checkpointOut": 10000, "synopsisOut": 1, "defectTracker": 0, "defectsOut": a.defout, "dirSOut": a.dirout, "flowOut": a.flowout}
os.makedirs(a.outdir, exist_ok=True); json.dump(inp, open(os.path.join(a.outdir, "input.json"), "w"), indent=4)
print(f"{a.outdir}: A {A:.2f} B {B:.2f} tip radius {B*B/A:.1f} domain {Lx}x{Ly} pop {inp['species'][0]['pop']}")
