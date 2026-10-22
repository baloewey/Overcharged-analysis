#!/usr/bin/env python3
"""Layer stationarity vs crowding (user, 2026-10-22: 'elongated layers were more diffuse, defects less stationary (Fourier); does crowding make them more
stationary?'). runs/cavity/wall_fourier.analyse on all 675 shape-sweep runs (band 4 cells, second half). At the lattice mode m = 2k-2: static fraction
f = |<rho_m>|^2/<|rho_m|^2> (1 = pinned lattice, 0 = sliding/disordered), S(m) = <|rho_m|^2>/<N>, and the decay time of the fluctuating part. -> shape_fourier.csv"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); import wall_fourier as WF  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        r = WF.analyse(run, 4.0, 0.5, mmax=60); k = int(r["k"]); i = 2 * k - 3
        x = None; return dict(run=os.path.relpath(run, SW), k=k, Nwall=r["Nwall"], f_lat=r["fstat"][i], S_lat=r["S"][i], mstar=r["mstar"], fstar=r["fstar"], Sstar=r["Sstar"], tau=r["tau"])
    except Exception as ex:
        return dict(run=os.path.relpath(run, SW), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "k", "Nwall", "f_lat", "S_lat", "mstar", "fstar", "Sstar", "tau", "err"]
    with open(os.path.join(H, "shape_fourier.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("errors", sum(1 for r in res if "err" in r))
