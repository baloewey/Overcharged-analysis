#!/usr/bin/env python3
"""Decomposition of the topological current across winding channels (user, 2026-10-14: is the bulk flux carried by defects created in the bulk rather than by
defects crossing the channel?). Charge conservation with pair processes: J_total = J_motion + J_events, where
  J_motion = (1/L) < sum_i q_i v_y,i >  over tracked defects in the bulk (> 4 cells from both walls), split by sign and by ORIGIN of the defect
             (born in a bulk pair / left a wall band (ejected) / other), and
  J_events = (1/L) sum over pair creations of q_+ (y_+ - y_-) + q_- ... = (1/L) sum_creations 0.5 (y_+ - y_-)  minus the same for annihilations, per unit time
             (a pair born with its +1/2 above its -1/2 carries charge upwards at birth; an annihilation of a +1/2 above a -1/2 carries it back down).
Also the mid-line crossings J_cross (linked defects only). All per unit length per 1000 t.u.; +y = from the bottom wall to the top wall. Second half of each run.
Tracking: kickout.kick_identity.track (same-sign links, gaps <= 3 frames); the x direction is periodic (links across the seam are missed; small).
Usage: current_decomp.py RUNDIR [...] -> stdout, current_decomp.txt"""
import os, sys, json, numpy as np
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "kickout"))
from edge_criteria import read_with_angle  # noqa
from kick_identity import track  # noqa


def one(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    cb, ct = J["BC"][2].get("customAnch", 0), J["BC"][3].get("customAnch", 0); ym = 0.5 * (y0 + y1)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); t = np.array(t); dt = float(np.median(np.diff(t))); T = len(fr) * dt
    D = [np.minimum(f[:, 1] - y0, y1 - f[:, 1]) if len(f) else np.zeros(0) for f in fr]
    ids, info = track(fr, D, t)
    origin = {}                                              # id -> 'pair' / 'wall' / 'other'
    for k, v in info.items(): origin[k] = {"PAIR": "pair", "WALL": "wall"}.get(v["type"], "other")
    was_wall = {}
    Jm = {(s, o): 0.0 for s in (1, -1) for o in ("pair", "wall", "other")}; cross = 0.0; ev_c = 0.0; ev_a = 0.0; nc = na = 0
    prev = {}
    for n, (f, idn) in enumerate(zip(fr, ids)):
        cur = {}
        for i, k in enumerate(idn):
            cur[k] = (f[i, 0], f[i, 1], np.sign(f[i, 2]), D[n][i])
            if D[n][i] < 4: was_wall[k] = True
        for k, (x, y, s, d) in cur.items():
            if k in prev:
                xp, yp, sp, dp = prev[k]
                if abs(x - xp) > L / 2: continue
                if d >= 4 and dp >= 4:
                    o = "wall" if was_wall.get(k) and origin[k] != "pair" else origin[k]
                    Jm[(int(s), o)] += 0.5 * s * (y - yp)
                if yp < ym <= y: cross += 0.5 * s
                elif y < ym <= yp: cross -= 0.5 * s
        # events: births and deaths of +/- pairs within 3 cells in the bulk
        born = [k for k in cur if k not in prev and n > 0]; died = [k for k in prev if k not in cur]
        def pairs(ks, pos):
            used, out = set(), []
            for a in ks:
                if pos[a][2] <= 0 or a in used: continue
                best = None
                for b in ks:
                    if b in used or pos[b][2] >= 0: continue
                    dx = (pos[b][0] - pos[a][0] + L / 2) % L - L / 2; dd = np.hypot(dx, pos[b][1] - pos[a][1])
                    if dd < 3 and (best is None or dd < best[0]): best = (dd, b)
                if best: used |= {a, best[1]}; out.append((pos[a][1], pos[best[1]][1]))
            return out
        for yp_, ym_ in pairs(born, cur): ev_c += 0.5 * (yp_ - ym_); nc += 1
        for yp_, ym_ in pairs(died, prev): ev_a -= 0.5 * (yp_ - ym_); na += 1
        prev = cur
    k = 1000 / (L * T)
    return dict(run=os.path.relpath(run, R_), cb=cb, ct=-ct if ct else 0, Jm={key: v * k for key, v in Jm.items()}, cross=cross * k, ev_c=ev_c * k, ev_a=ev_a * k, nc=nc, na=na)


if __name__ == "__main__":
    runs = [os.path.abspath(a) for a in sys.argv[1:] if os.path.isdir(a)]
    from multiprocessing import Pool
    with Pool(min(6, len(runs))) as pool: R = pool.map(one, runs)
    L = ["TOPOLOGICAL CURRENT decomposition (per unit length per 1000 t.u.; +y = bottom wall -> top wall). Wall charges (bottom/top) from customAnch.",
         "run | walls | motion: +1/2 pair-born, +1/2 ejected, +1/2 other | -1/2 pair-born, -1/2 other | motion total | events: creations, annihilations | TOTAL (motion+events) | mid-line crossings"]
    for r in R:
        m = r["Jm"]; mt = sum(m.values()); ev = r["ev_c"] + r["ev_a"]
        L.append(f"{r['run']:46s} | {r['cb']:+d}/{-r['ct'] if r['ct'] else 0:+d} | {m[(1,'pair')]:+.3f} {m[(1,'wall')]:+.3f} {m[(1,'other')]:+.3f} | {m[(-1,'pair')]:+.3f} {m[(-1,'other')]+m[(-1,'wall')]:+.3f} | "
                 f"{mt:+.3f} | {r['ev_c']:+.3f} {r['ev_a']:+.3f} | {mt + ev:+.3f} | {r['cross']:+.3f}")
    print("\n".join(L)); open(os.path.join(H, os.environ.get("OUT", "current_decomp.txt")), "w").write("\n".join(L) + "\n")
