#!/usr/bin/env python3
"""
Wall +1/2 travel along the wall (master document, Sec. on the wall lattice in ovals and ellipses).
Gap-closed +1/2 tracks (build_tracks, rlink 8, gap 4), frames after 20%. For every pair of frames lag apart (50 and 500 time
units) in which a track stays within 4 cells of the wall throughout, the unwrapped along-wall displacement |ds| is recorded.
Reported: number of such pairs, median |ds|, fraction with |ds| > half the site spacing P/(2(2k-2)).
Run from runs/:  python3 ../reports/master/scripts/wall_travel.py > ../reports/master/data/wall_travel.txt
"""
import sys, json, numpy as np
import sys, json, numpy as np
from multiprocessing import Pool
sys.path.insert(0,"cavity")
from edge_criteria import read_with_angle
from wall_fourier import wall_polyline, project
from defect_tracks import build_tracks
BAND=4.0
def one(arg):
    lab,r=arg
    bc=json.load(open(r+"/input.json"))["BC"][0]; k=int(bc["customAnch"])
    W,s,P=wall_polyline(bc)
    t,fr=read_with_angle(r,np.array(bc["Q"][:2],float),0.2); dt=t[1]-t[0]
    tr,_=build_tracks(fr,dt,8.0,4)
    spacing=P/max(2*k-2,1)
    res={}
    stay_frac=[]
    for lag_t in (50,500):
        lag=int(round(lag_t/dt)); ds=[]
        for T in tr.values():
            if T["q"]<0 or len(T["n"])<lag+1: continue
            d,sa=project(W,s,T["xy"])
            n=T["n"]; inw=d<BAND
            su=np.unwrap(sa/P*2*np.pi)*P/(2*np.pi)       # unwrapped arclength
            idx={v:i for i,v in enumerate(n)}
            for i,ni in enumerate(n):
                j=idx.get(ni+lag)
                if j is None: continue
                if inw[i:j+1].all(): ds.append(abs(su[j]-su[i]))
        ds=np.array(ds)
        res[lag_t]=(len(ds), np.median(ds) if len(ds) else np.nan, np.mean(ds>0.5*spacing) if len(ds) else np.nan)
    return lab,k,spacing,res
C="cavity/server_main/main_act0.07"
runs=[("circle 1.43 la s1",f"{C}/rho0.06/k06/s01"),("circle 1.43 la s2",f"{C}/rho0.06/k06/s02"),("circle 2.15 la s1",f"{C}/rho0.04/k06/s01"),
      ("passive ellipse 1.5","cassini/P5_passive_refs/ellipse_eps1.5_k06"),("passive ellipse 2.0","cassini/P5_passive_refs/ellipse_eps2.0_k06")]
for e in ("1.2","1.5","2.0"):
    for sd in (1,2): runs.append((f"ellipse {e} s{sd}",f"cassini/W3_ellipse/eps{e}_k06/s0{sd}"))
for b in ("2.5","1.6","1.3"): runs.append((f"oval H34 b/a {b} s1",f"cassini/P3/H34_ba{b}_k06/s01"))
with Pool(12) as p:
    out=p.map(one,runs)
print(f"{'run':22s} {'k':>2} {'half spacing':>12} | lag 50: n  median|ds|  frac>half-spacing | lag 500: n  median|ds|  frac>half-spacing")
for lab,k,sp,res in out:
    a,b=res[50],res[500]
    print(f"{lab:22s} {k:2d} {sp/2:12.1f} | {a[0]:7d} {a[1]:8.2f} {a[2]:8.3f}        | {b[0]:7d} {b[1]:8.2f} {b[2]:8.3f}")
