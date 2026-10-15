#!/usr/bin/env python3
"""How are the wall layers of a mixed-winding channel replenished? (user, 2026-10-15). Charge budget of each wall band (< 4 cells from the wall), from tracked
identities (kickout.kick_identity.track; second half; rates per unit length per 1000 t.u.):
  * crossings of the band edge by +1/2 and by -1/2, in and out;
  * pair creations inside the band (and where each member goes: stays in the band or leaves for the bulk);
  * annihilations inside the band, classified by the members' origin (both band residents / one arrived from the bulk);
  * defects appearing or vanishing alone (tracking noise).
Steady state check: for each band, d(charge)/dt = 0.5*(+in - +out) - 0.5*(-in - -out) + 0 (pair events are neutral) ~ 0.
Sign convention of the current J (+y = from the bottom wall to the top wall): q v_y, i.e. +1/2 moving up and -1/2 moving down both count positive.
Usage: layer_budget.py RUNDIR [...] -> layer_budget.txt"""
import os, sys, json, numpy as np
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "kickout"))
from edge_criteria import read_with_angle  # noqa
from kick_identity import track  # noqa


def one(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); t = np.array(t); dt = float(np.median(np.diff(t))); T = len(fr) * dt
    reg = lambda y: 0 if y - y0 < 4 else (2 if y1 - y < 4 else 1)          # 0 bottom band, 1 bulk, 2 top band
    D = [np.minimum(f[:, 1] - y0, y1 - f[:, 1]) if len(f) else np.zeros(0) for f in fr]
    ids, info = track(fr, D, t)
    C = {}                                                                     # counters
    def add(k, v=1): C[k] = C.get(k, 0) + v
    prev = {}
    for n, (f, idn) in enumerate(zip(fr, ids)):
        cur = {k: (f[i, 0], f[i, 1], int(np.sign(f[i, 2])), reg(f[i, 1])) for i, k in enumerate(idn)}
        for k, (x, y, s, r) in cur.items():
            if k in prev and prev[k][3] != r:
                rp = prev[k][3]
                for band in (0, 2):
                    if rp != band and r == band: add((band, s, "in"))
                    if rp == band and r != band: add((band, s, "out"))
        if n > 0:
            born = [k for k in cur if k not in prev]; died = [k for k in prev if k not in cur]
            def pairs(ks, pos):
                used, out = set(), []
                for a in ks:
                    if pos[a][2] <= 0 or a in used: continue
                    for b in ks:
                        if b in used or pos[b][2] >= 0: continue
                        if np.hypot((pos[b][0] - pos[a][0] + L / 2) % L - L / 2, pos[b][1] - pos[a][1]) < 3: used |= {a, b}; out.append((a, b)); break
                return out, set(ks) - used
            pc, lone_b = pairs(born, cur); pa, lone_d = pairs(died, prev)
            for a, b in pc:
                r = reg(0.5 * (cur[a][1] + cur[b][1])); add((r, "create"))
            for a, b in pa:
                r = reg(0.5 * (prev[a][1] + prev[b][1])); add((r, "annihilate"))
            for k in lone_b: add((cur[k][3], cur[k][2], "appear"))
            for k in lone_d: add((prev[k][3], prev[k][2], "vanish"))
        prev = cur
    k = 1000 / (L * T); return {key: v * k for key, v in C.items()}, J


if __name__ == "__main__":
    runs = [os.path.abspath(a.rstrip("/")) for a in sys.argv[1:] if os.path.isdir(a)]
    from multiprocessing import Pool
    with Pool(min(8, len(runs))) as pool: R = pool.map(one, runs)
    L = ["WALL-BAND BUDGET (per unit length per 1000 t.u.; band = < 4 cells from a wall; in/out = crossings of the band edge).",
         "run | band (wall charge) | +1/2 in, out | -1/2 in, out | pair creations, annihilations in the band | lone appear/vanish (+, -) | net charge flux into the band by crossings"]
    for run, (C, J) in zip(runs, R):
        g = lambda *k: C.get(tuple(k), 0.0)
        for band, lab in ((0, f"bottom ({J['BC'][2].get('customAnch',0):+d})"), (2, f"top ({-J['BC'][3].get('customAnch',0):+d})")):
            net = 0.5 * (g(band, 1, "in") - g(band, 1, "out")) - 0.5 * (g(band, -1, "in") - g(band, -1, "out"))
            L.append(f"{os.path.relpath(run, R_):42s} | {lab:12s} | {g(band,1,'in'):.3f} {g(band,1,'out'):.3f} | {g(band,-1,'in'):.3f} {g(band,-1,'out'):.3f} | "
                     f"{g(band,'create'):.3f} {g(band,'annihilate'):.3f} | {g(band,1,'appear'):.3f}/{g(band,1,'vanish'):.3f} {g(band,-1,'appear'):.3f}/{g(band,-1,'vanish'):.3f} | {net:+.4f}")
        L.append(f"{'':42s} | bulk         | pair creations {C.get((1,'create'),0):.3f}, annihilations {C.get((1,'annihilate'),0):.3f}")
    print("\n".join(L)); open(os.path.join(H, os.environ.get("OUT", "layer_budget.txt")), "w").write("\n".join(L) + "\n")
