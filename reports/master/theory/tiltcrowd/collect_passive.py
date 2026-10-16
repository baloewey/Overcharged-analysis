#!/usr/bin/env python3
"""Collect per-run flow/depth/tilt data of the passive imposed-vortex runs (runs/swirl_crowd; swirl 0.01 'k*' and 0.03 'A03_k*', k <= 13 where the ring ejects)
with the same functions as the active set (defect_flow.one, wallshear.one). -> passive_runs.npz (for the tilt report)."""
import os, sys, glob, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import defect_flow as DF, wallshear as WS
def one(run):
    a = DF.one(run); w = WS.one(run); c = a["core"]; p = a["prof"][2.0]
    return dict(name=os.path.basename(run), A=0.03 if os.path.basename(run).startswith("A03") else 0.01, k=a["k"], a=a["a"], d=a["depth"], tilt=a["tilt"],
                s=w["e_band"] + w["w_band"], vth=w["vth"], phi_core=c["phi"], exx_core=c["exx"], ux_core=c["ux"], phi2=p["phi"], tau2=p["tau_pred"])
if __name__ == "__main__":
    runs = sorted(r for r in glob.glob("/home/baloewe/Physics/Overcharged/runs/swirl_crowd/*k*_s*") if os.path.exists(r + "/DONE"))
    runs = [r for r in runs if int(os.path.basename(r).split("k")[1][:2]) <= (13 if "A03" in r else 8)]
    with Pool(6) as p: R = p.map(one, runs)
    np.savez(os.path.join(H, "passive_runs.npz"), **{key: np.array([x[key] for x in R]) for key in R[0]}); print(len(R), "runs")
