#!/usr/bin/env python3
"""Continuum core rotation where the defect sits deeper (2026-10-22): flat wall, a_s 30, 40, 60 (deeper cores), s = 0 (control), 0.01, 0.03.
The estimator is biased when the core overlaps the wall (core_rotation_controls: -9 to -23 deg at rest for depth 0.7-1.4); compare sheared with the
rest control at similar depth. -> core_rotation_deep.txt"""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore"); sys.path.insert(0, "."); import core_rotation as CR
from multiprocessing import Pool
if __name__ == "__main__":
    jobs = [(a_s, s, 2.0) for a_s in (30.0, 40.0, 60.0) for s in (0.0, 0.01, 0.03)]
    with Pool(9) as p: C = p.map(CR.continuum, jobs)
    L = ["CONTINUUM flat wall, deeper cores: a_s | s | tilt | depth | ring 2-3 | ring 3-4.5 (deg; s = 0 is the no-flow control)"]
    for x in C: L.append(f"{x['job'][0]:4.0f} | {x['job'][1]:.2f} | {x['tilt']:6.1f} | {x['depth']:.2f} | {x['r23']:+6.1f} | {x['r345']:+6.1f}")
    open("core_rotation_deep.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
