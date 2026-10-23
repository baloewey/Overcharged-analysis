#!/usr/bin/env python3
"""Channel crossover in wall density (user, 2026-10-22: 'run channels at higher q; the leak should switch off near rho 0.08-0.1'). W 40, L 192, both walls
hosting +1/2 layers (eMode-4 winding bottom +m, top -m with the normal flipped): runs/channel_qcross (m 6, 7, 8, 10, 11; 2 seeds), channel_qscan (m 3, 9, 12),
channel/W40_m6_L192, channel_fix crowd sym (m 18, 24; tracker blind spot, unreliable). Second half of each run. rho = 2m/L per wall.
With the validated tracker (NEW v2, crowd_protect/tracker_validation/trackers.track_new, gate 9, periodic x):
  p_esc (shape_runs definition: wall +1/2 at d < 4 every 10 t.u., deeper than 6 cells 20 t.u. later); clean escape hazard per wall defect and entry hazard
  per bulk +1/2 (exchange definitions, per 1000 t.u.); tracking-free: wall +1/2 per wall / 2m (occupancy), bulk charge per length, and the mid-line
  current = slope of Q_up(t) (must be ~0). -> channel_crossover.txt"""
import glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../../crowd_protect/tracker_validation"); sys.path.insert(0, "../../crowd_protect/exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs"; NW = 10
def one(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1); m = abs(J["BC"][2]["customAnch"])
    t, fr = X.read_with_angle(run, np.array([0.0, 0.0]), 0.5); nf = len(fr); dt = t[1] - t[0]; dist = lambda xy: np.minimum(xy[:, 1] - y0, y1 - xy[:, 1])
    T = track_new(fr, 9.0, L=L); st = es = 0; esc = ent = 0; expW = expB = 0.0
    for tr in T:
        if tr["sg"] < 0: continue
        nn = np.array(tr["n"]); d = dist(np.array(tr["xy"])); pos = dict(zip(nn, d))
        for n, dv in pos.items():
            if n % 5 or n + 10 >= nf or dv >= 4: continue
            st += 1; e = pos.get(n + 10); es += e is not None and e > 6
        full = np.full(nn[-1] - nn[0] + 1, np.nan); full[nn - nn[0]] = d
        for i in range(1, len(full)):
            if np.isnan(full[i]): full[i] = full[i - 1]
        s = np.where(full < 4, 1, np.where(full > 6, 2, 0)); runW = runB = 0
        for a in range(len(s)):
            if s[a] == 1:
                if runB >= NW and a + NW < len(s) and np.all(s[a:a + NW + 1] == 1): ent += 1
                runB = 0; runW += 1; expW += dt if runW > NW else 0
            elif s[a] == 2:
                if runW >= NW and a + NW < len(s) and not np.any(s[a:a + NW + 1] == 1): esc += 1
                runW = 0; runB += 1; expB += dt if runB > NW else 0
    nwall = np.mean([((dist(f[:, :2]) < 4) & (f[:, 2] > 0)).sum() / 2 for f in fr if len(f)])
    qb = np.mean([0.5 * np.sign(f[dist(f[:, :2]) > 4, 2]).sum() for f in fr if len(f)])
    qu = np.array([0.5 * np.sign(f[f[:, 1] > ym, 2]).sum() for f in fr if len(f)]); J_ = 1000 * np.polyfit(np.arange(len(qu)) * dt, qu, 1)[0] / L
    return dict(run=os.path.relpath(run, R_), m=m, rho=2 * m / L, p_esc=es / max(st, 1), h_esc=1000 * esc / max(expW, 1), h_ent=1000 * ent / max(expB, 1),
                occ=nwall / (2 * m), qb_len=qb / (2 * L), J=J_)
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(R_ + "/channel_qcross/W40_*/DONE") + glob.glob(R_ + "/channel_qscan/W40_*/DONE")
                  + glob.glob(R_ + "/channel/W40_m6_L192_s*/DONE") + glob.glob(R_ + "/channel_fix/crowd_W40_*_sym_s*/DONE"))
    with Pool(8) as p: res = p.map(one, runs)
    by = {}
    for r in res: by.setdefault(r["m"], []).append(r)
    Lt = ["CHANNEL CROSSOVER (W 40, L 192; validated tracker; means over the runs at each m, all seeds/sets pooled)",
          "m | rho = 2m/L | runs | p_esc per wall +1/2 per 20 t.u. | clean escape hazard /1000 t.u. | entry hazard /1000 t.u. | wall occupancy N_wall/(2m) | bulk charge per wall length | mid-line current (must be ~0)"]
    for m in sorted(by):
        a = by[m]; mn = lambda k: np.mean([x[k] for x in a])
        Lt.append(f"{m:2d} | {2 * m / 192:.3f} | {len(a)} | {mn('p_esc'):.4f} | {mn('h_esc'):6.2f} | {mn('h_ent'):5.2f} | {mn('occ'):.3f} | {mn('qb_len'):.4f} | {mn('J'):+.4f}" + ("   (m >= 18: tracker blind spot, unreliable)" if m >= 18 else ""))
    open("channel_crossover.txt", "w").write("\n".join(Lt) + "\n"); print("\n".join(Lt))
