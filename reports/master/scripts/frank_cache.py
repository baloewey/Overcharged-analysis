#!/usr/bin/env python3
"""Run frank_constant.analyse once on every Frank-calibration box and cache the results for the figure scripts (2026-10-21).
runs/Kcalib (production coupling: tumble 2, shearSusc 0.5; L64 s1/s2, L96 s1) and runs/Kcalib_lambda (L64, 2 seeds each: shearSusc 0; tumble 0, 1, -2).
-> data/frank_cache.pkl  {run name: analyse() dict}"""
import sys, pickle
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/Kcalib"); import frank_constant as FC
R = "/home/baloewe/Physics/Overcharged/runs/"
RUNS = [f"Kcalib/{r}" for r in ("L64_s1", "L64_s2", "L96_s1")] + [f"Kcalib_lambda/L64_{c}_s{s}" for c in ("susc0", "tum0", "tum1", "tumm2") for s in (1, 2)]
def one(r): return r, FC.analyse(R + r, 0.8, 200.0)
if __name__ == "__main__":
    with Pool(6) as p: res = dict(p.map(one, RUNS))
    pickle.dump(res, open("data/frank_cache.pkl", "wb")); print({k: (round(v["K1"], 2), round(v["K3"], 2)) for k, v in res.items()})
