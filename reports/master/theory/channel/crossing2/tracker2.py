#!/usr/bin/env python3
"""Birth-to-death defect trajectories in a periodic channel and their mid-line crossings (user, 2026-10-22).
Per sign, frame to frame (2 t.u.): predicted position = last position + last displacement (zero for a new track); global assignment (Hungarian,
scipy linear_sum_assignment) on the periodic distance between predictions and detections, gated at GATE cells; unmatched tracks survive up to GAP missed
frames (prediction carried on); unmatched detections start new tracks. A track's crossings: every sign change of y - y_mid between consecutive positions
(gaps included). Pairs: a +1/2 and a -1/2 born (ending) in the same frame within PAIR cells on opposite sides of the mid-line move charge across it and are
counted as such. Check: cumulative (crossings + straddling pairs) against the tracking-free Q_up(t) - Q_up(0).
API: build(run) -> dict(tracks, t, ym, L, Qup, cum, parts)."""
import json, sys, numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
GATE, GAP, PAIR = 6.0, 2, 5.0
def build(run, gate=GATE, gap=GAP, predict=True):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); nf = len(fr)
    tracks = []; active = {+1: [], -1: []}                      # active: list of track indices
    for n, f in enumerate(fr):
        for sg in (+1, -1):
            p = f[np.sign(f[:, 2]) == sg][:, :2] if len(f) else np.zeros((0, 2)); act = active[sg]; used = np.zeros(len(p), bool); keep = []
            if act and len(p):
                pred = []
                for ti in act:
                    tr = tracks[ti]; x, y = tr["xy"][-1]; m = n - tr["n"][-1]
                    v = (np.array(tr["xy"][-1]) - np.array(tr["xy"][-2])) / (tr["n"][-1] - tr["n"][-2]) if predict and len(tr["n"]) > 1 else np.zeros(2)
                    if predict and len(tr["n"]) > 1: v[0] = ((v[0] + L / 2) % L) - L / 2
                    pred.append((x + v[0] * m, y + v[1] * m))
                pred = np.array(pred); dx = (p[None, :, 0] - pred[:, None, 0] + L / 2) % L - L / 2; dy = p[None, :, 1] - pred[:, None, 1]; C = np.hypot(dx, dy)
                r, c = linear_sum_assignment(np.where(C < gate, C, 1e6))
                for a, b in zip(r, c):
                    if C[a, b] < gate:
                        tr = tracks[act[a]]; tr["n"].append(n); tr["xy"].append(tuple(p[b])); used[b] = True
            for ti in act:
                if n - tracks[ti]["n"][-1] <= gap: keep.append(ti)
            for b in np.flatnonzero(~used):
                tracks.append(dict(sg=sg, n=[n], xy=[tuple(p[b])])); keep.append(len(tracks) - 1)
            active[sg] = keep
    # crossings along tracks
    cum = np.zeros(nf); parts = dict(p_up=0, p_down=0, m_up=0, m_down=0, pair_birth=0.0, pair_death=0.0)
    for tr in tracks:
        y = np.array(tr["xy"])[:, 1]; s = np.sign(y - ym); s[s == 0] = 1; q = 0.5 * tr["sg"]
        for i in np.flatnonzero(s[1:] != s[:-1]):
            up = s[i + 1] > 0; cum[tr["n"][i + 1]] += q if up else -q
            parts[("p_" if tr["sg"] > 0 else "m_") + ("up" if up else "down")] += 1
        tr["cross"] = int(np.sum(s[1:] != s[:-1]))
    # pairs straddling the mid-line
    for key, idx in (("pair_birth", 0), ("pair_death", -1)):
        ev = {}
        for tr in tracks:
            nn = tr["n"][idx]
            if (idx == 0 and nn == 0) or (idx == -1 and nn >= nf - 1 - gap): continue
            ev.setdefault(nn, []).append((tr["sg"], tr["xy"][idx]))
        for nn, es in ev.items():
            P = [e for e in es if e[0] > 0]; M = [e for e in es if e[0] < 0]; usedm = set()
            for sg, (x, y) in P:
                best = None
                for j, (_, (x2, y2)) in enumerate(M):
                    if j in usedm: continue
                    d = np.hypot(((x2 - x + L / 2) % L) - L / 2, y2 - y)
                    if d < PAIR and (best is None or d < best[0]): best = (d, j, y2)
                if best:
                    usedm.add(best[1]); y2 = best[2]
                    if (y > ym) != (y2 > ym):
                        dq = (0.5 if y > ym else -0.5) * (1 if idx == 0 else -1)          # birth: + above adds +1/2 above ... death removes it
                        cum[nn if idx == 0 else min(nn + 1, nf - 1)] += dq; parts[key] += dq
    qup = np.array([0.5 * np.sign(f[f[:, 1] > ym, 2]).sum() if len(f) else np.nan for f in fr])
    return dict(tracks=tracks, t=t, ym=ym, L=L, y0=y0, y1=y1, Qup=qup, cum=np.cumsum(cum), parts=parts, J=J, fr=fr)
if __name__ == "__main__":
    run = sys.argv[1]
    for pr, g in ((False, 2.5), (True, 6.0), (True, 9.0)):
        B = build(run, gate=g, predict=pr); q = B["Qup"] - np.nanmean(B["Qup"][:25]); c = B["cum"]
        lens = np.array([len(tr["n"]) for tr in B["tracks"]])
        print(f"predict {pr}, gate {g}: tracks {len(B['tracks'])}, median length {np.median(lens) * 2:.0f} t.u.; crossings {B['parts']}; "
              f"tracked cumulative at end {c[-1]:+.1f} vs Q_up change {np.nanmean(q[-25:]):+.1f}; rms(tracked - Q_up) {np.sqrt(np.nanmean((c - q) ** 2)):.1f}", flush=True)
