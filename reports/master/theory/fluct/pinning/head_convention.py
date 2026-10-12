#!/usr/bin/env python3
"""Which way does a +1/2 self-propel relative to the recorded defect angle (direction of div Q = tail for a +1/2)?
Bulk +1/2 (more than 6 cells from the wall) in active circles (k = 6, act 0.07, three densities, seed 1): nearest-neighbour matching between
consecutive defect frames (dt = 2 t.u., match radius 3 cells, unambiguous matches only). Displacement d per frame; components along the
recorded angle e (cos psi, sin psi). Extensile: self-propulsion along the head = -e -> <d . e> < 0. Advection by the mean flow cancels on
average (orientation independent of the local flow direction only approximately; checked by the passive control: <d . e> ~ 0).
Writes head_convention.txt."""
import glob, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
for lab, run in [("active rho0.05 k06 s01", ROOT + "/cavity/server_main/main_act0.07/rho0.05/k06/s01"),
                 ("active rho0.06 k06 s01", ROOT + "/cavity/server_main/main_act0.07/rho0.06/k06/s01"),
                 ("active rho0.08 k06 s01", ROOT + "/cavity/server_main/main_act0.07/rho0.08/k06/s01"),
                 ("active ellipse eps1.5 s01", ROOT + "/cassini/W3_ellipse/eps1.5_k06/s01"),
                 ("passive rho0.06 k06 (control)", ROOT + "/cavity/server_main/main_passive/rho0.06/k06/s01")]:
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; a = bc.get("aInv", [1, 1])
    A, B = R * a[0], R * a[1]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.2)
    par, perp, wpar = [], [], []
    for f0, f1 in zip(fr[:-1], fr[1:]):
        p0 = f0[f0[:, 2] > 0]; p1 = f1[f1[:, 2] > 0]
        if len(p0) == 0 or len(p1) == 0: continue
        tr = cKDTree(p1[:, :2]); dist, j = tr.query(p0[:, :2])
        for k in range(len(p0)):
            if dist[k] > 3: continue
            if np.sum(np.hypot(*(p1[:, :2] - p0[k, :2]).T) < 3) > 1: continue          # ambiguous
            x, y = p0[k, :2]; inside = (x / A) ** 2 + (y / B) ** 2
            dvec = p1[j[k], :2] - p0[k, :2]; e = np.array([np.cos(p0[k, 3]), np.sin(p0[k, 3])])
            if np.sqrt(inside) < 1 - 6 / min(A, B): par.append(dvec @ e); perp.append(dvec @ np.array([-e[1], e[0]]))
            elif np.sqrt(inside) > 1 - 4 / min(A, B): wpar.append(dvec @ e)
    par, perp = np.array(par), np.array(perp)
    say(f"{lab:30s}: bulk +1/2 steps {len(par)}: <d.e> = {par.mean():+.3f} +- {par.std() / np.sqrt(len(par)):.3f} cells per 2 t.u. "
        f"(perp {perp.mean():+.3f}); fraction moving against e {np.mean(par < 0):.2f}; wall +1/2 <d.e> {np.mean(wpar) if wpar else np.nan:+.3f}")
open(os.path.join(H, "head_convention.txt"), "w").write("\n".join(L) + "\n")
