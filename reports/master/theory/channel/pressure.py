#!/usr/bin/env python3
"""Pressure-driven winding channels (user, 2026-10-16; runs/channel_pressure): body force g_x = +-0.001, +-0.003 on every particle, L 192, W 24 (check input),
mixed +6/-6 and symmetric +6/+6 active (2 seeds), passive mixed and homeotropic controls. Per run (second half): net flow <v_x> and its sign changes
(net_flow.one), wall-layer counts and bulk charge (crowd_an.one), and for mixed channels the topological current (crossings). Rectification test: does the
response to +g and -g differ beyond the sign (|<v_x>(+g)| vs |<v_x>(-g)|), and do the layers or the current depend on the flow direction? -> pressure.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import net_flow as NF, crowd_an as CA
R_ = "/home/baloewe/Physics/Overcharged/runs/channel_pressure"
def one(run):
    J = json.load(open(run + "/input.json")); g = (J.get("grav") or [0])[0]
    a = NF.one(run); b = CA.one(run)
    return dict(run=os.path.basename(run), g=g, W=b["W"], mixed=b["mixed"], act=J["species"][0]["act"], vx=a["mean"] if a else np.nan, vrms=a["vrms"] if a else np.nan, sc=a["sc"] if a else np.nan,
                fpos=a["fpos"] if a else np.nan, bp=b["bp"], bm=b["bm"], tp=b["tp"], tm=b["tm"], qb=b["qb"], cross=b.get("cross", np.nan), tilt=b["tilt"])
if __name__ == "__main__":
    runs = sorted(r.rstrip("/") for r in glob.glob(R_ + "/*/") if os.path.exists(r + "DONE"))
    with Pool(6) as p: res = p.map(one, runs)
    L = ["PRESSURE-DRIVEN WINDING CHANNELS (second half). g = body force along x; <v_x> channel-averaged net flow; layers within 4 cells; current = net mid-line crossings per length per 1000 t.u.",
         "run | g | W | act | <v_x> | frac>0 | sign changes/1000 t.u. | v_rms | bottom +1/2 -1/2 | top +1/2 -1/2 | bulk charge | current | tilt bottom layer"]
    for r in sorted(res, key=lambda r: (r["run"].split("_g")[0], r["act"], abs(r["g"]), r["g"], r["run"])):
        L.append(f"{r['run']:24s} | {r['g']:+.3f} | {r['W']:.0f} | {r['act']:.2f} | {r['vx']:+.3f} | {r['fpos']:.2f} | {r['sc']:5.1f} | {r['vrms']:.3f} | {r['bp']:5.2f} {r['bm']:4.2f} | {r['tp']:4.2f} {r['tm']:5.2f} | {r['qb']:+.2f} | "
                 + (f"{r['cross']:+.4f}" if r["mixed"] else "-") + f" | {r['tilt']:5.1f}")
    open(os.path.join(H, "pressure.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
