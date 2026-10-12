#!/usr/bin/env python3
"""
Check how the imposed winding is distributed along an elliptical wall for eMode 0/1/2.

On the contour (a-d) cos t, (b-d) sin t just inside the wall, sample the director (Q-tensor bilinear
interpolation from directorfield.dat), accumulate the nematic winding, and report
  - total winding (should be k)
  - fraction of the winding that happens in the tip regions |x| > a/sqrt(2)
  - predictions for that fraction: normal-angle (eMode 0), parameter-uniform (1), arclength-uniform (2)
Also runs the charge audit on defects.dat.
Usage: check_emode.py RUN [RUN ...] [--d 1.2]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cavityValidation"))
sys.path.insert(0, os.path.join(HERE, "..", "cavity"))
from cavity_check import bilinear, read_director  # noqa: E402
from defect_audit import audit_frames, audit_summary, read_frames  # noqa: E402


def predictions(a, b, tip):
    t = np.linspace(0, 2 * np.pi, 200001)[:-1]
    x = a * np.cos(t)
    ds = np.sqrt((a * np.sin(t))**2 + (b * np.cos(t))**2)
    th_n = np.unwrap(np.arctan2(a * np.sin(t), b * np.cos(t)))      # outward normal angle
    dth = np.gradient(th_n)
    m = np.abs(x) > tip
    return {0: dth[m].sum() / dth.sum(), 1: m.mean(), 2: ds[m].sum() / ds.sum()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--d", type=float, default=1.2)
    a_ = ap.parse_args()
    print(f"{'run':<14} {'eMode':>5} {'total winding':>13} {'tip fraction':>12} | predicted e0 / e1 / e2 | audit q=k  N")
    for run in a_.runs:
        inp = json.load(open(os.path.join(run, "input.json")))
        bc = inp["BC"][0]
        L = tuple(inp["domain"][:2]); C = bc["Q"][:2]
        a, b = bc["R"] * bc["aInv"][0], bc["R"] * bc["aInv"][1]
        k, e = bc["customAnch"], bc.get("eMode", 0)
        t_, Qxx, Qxy, S = read_director(os.path.join(run, "directorfield.dat"), L)
        keep = t_ >= 0.5 * t_.max()
        tt = np.linspace(0, 2 * np.pi, 1441)[:-1]
        xs, ys = C[0] + (a - a_.d) * np.cos(tt), C[1] + (b - a_.d) * np.sin(tt)
        tipm = np.abs((a - a_.d) * np.cos(tt)) > (a - a_.d) / np.sqrt(2)
        tot, tipf = [], []
        for f in np.nonzero(keep)[0]:
            psi = 0.5 * np.arctan2(bilinear(Qxy[f], xs, ys), bilinear(Qxx[f], xs, ys))
            dpsi = np.diff(np.r_[psi, psi[0]])
            dpsi = (dpsi + np.pi / 2) % np.pi - np.pi / 2
            tot.append(dpsi.sum() / (2 * np.pi))
            tipf.append(dpsi[tipm].sum() / dpsi.sum() if abs(dpsi.sum()) > 1e-9 else np.nan)
        pr = predictions(a - a_.d, b - a_.d, (a - a_.d) / np.sqrt(2))
        _, fr = read_frames(os.path.join(run, "defects.dat"))
        au = audit_summary(audit_frames(fr[len(fr) // 2:], k))
        N = np.mean([len(x) for x in fr[len(fr) // 2:]])
        print(f"{os.path.basename(run.rstrip('/')):<14} {e:>5} {np.mean(tot):13.2f} {np.nanmean(tipf):12.3f} | "
              f"{pr[0]:.3f} / {pr[1]:.3f} / {pr[2]:.3f}   | {au['frac_ok']:6.2f} {N:5.1f}")


if __name__ == "__main__":
    main()
