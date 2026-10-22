#!/usr/bin/env python3
"""Synthetic ground-truth test of the three linkers (2026-10-22). Circle R 31.83 (P 200). Wall layer: N_w +1/2 bound to evenly spaced sites at depth
1.2-1.8 with per-frame jitter (sigma 0.6 cells along, 0.3 across). Bulk +1/2: step lengths drawn from the measured bulk step distribution
(steps_bulk_*.npy, new tracker on real runs), direction persistent (angle noise 0.5 rad/frame), reflected at d < 0.5. Wall defects relax towards the midpoint of their neighbours (factor 0.2 per frame; the layer repels). Escapes: a wall defect leaves at
rate r per frame, moving inward (persistence 0.15 rad) and becomes a bulk defect; a bulk defect reaching d < 3 joins the wall in the middle of the nearest gap (entry).
Detection: positions rounded to 0.5 cells (cluster means), 2% of detections missed per frame. 2 t.u. per frame, 2500 frames.
Truth vs each linker (NEW gate 9; NN 3.5; CHAIN 10 t.u. / 6 cells):
  p_esc (shape_runs definition: start d < 4 every 5 frames, deeper than 6 cells 10 frames later); clean escapes and entries (exchange definitions);
  identity switches per 1000 links (bulk and wall separately); track fragmentation (tracks per true bulk trajectory).
Scenarios: sparse (N_w 8, r 0.015) and crowded (N_w 28, r 0.0012), bulk population ~5; 3 seeds. -> synth.txt"""
import numpy as np
from trackers import track_new, track_nn, chain_follow
R = 31.830989; P = 2 * np.pi * R
def simulate(Nw, r, seed, nf=2500, steps=None):
    rng = np.random.default_rng(seed); sites = 2 * np.pi * np.arange(Nw) / Nw
    P_ = []                                                           # particles: dict(state 'w'/'b', theta/site, d, pos, ang)
    for i in range(Nw): P_.append(dict(id=i, st="w", phi=sites[i], d=1.5))
    for i in range(5):
        a = rng.uniform(0, 2 * np.pi); rr = rng.uniform(0, R - 8); P_.append(dict(id=Nw + i, st="b", x=rr * np.cos(a), y=rr * np.sin(a), ang=rng.uniform(0, 2 * np.pi)))
    nid = Nw + 5; truth = []; frames = []
    for n in range(nf):
        rows = []; tr = []
        W_ = sorted([p for p in P_ if p["st"] == "w"], key=lambda p: p["phi"] % (2 * np.pi))     # elastic relaxation of the wall layer
        if len(W_) > 2:
            ph = np.array([p["phi"] % (2 * np.pi) for p in W_]); prv = np.roll(ph, 1); nxt = np.roll(ph, -1); prv[0] -= 2 * np.pi; nxt[-1] += 2 * np.pi
            for p, a, b, c in zip(W_, prv, ph, nxt): p["phi"] = b + 0.2 * (0.5 * (a + c) - b)
        for p in P_:
            if p["st"] == "w":
                if rng.random() < r:
                    p["st"] = "b"; x0, y0 = (R - p["d"]) * np.cos(p["phi"]), (R - p["d"]) * np.sin(p["phi"]); p.update(x=x0, y=y0, ang=p["phi"] + np.pi, esc=12)
                else:
                    ph = p["phi"] + rng.normal(0, 0.6) / R; d = np.clip(1.5 + rng.normal(0, 0.3), 0.8, 2.6); x, y = (R - d) * np.cos(ph), (R - d) * np.sin(ph)
            if p["st"] == "b":
                s = rng.choice(steps); pers = 0.15 if p.get("esc", 0) > 0 else 0.5; p["esc"] = p.get("esc", 0) - 1
                p["ang"] += rng.normal(0, pers); x, y = p["x"] + s * np.cos(p["ang"]), p["y"] + s * np.sin(p["ang"]); rr = np.hypot(x, y)
                if rr > R - 0.5:                                      # reflect
                    x, y = x * (R - 0.5) / rr, y * (R - 0.5) / rr; p["ang"] = np.arctan2(y, x) + np.pi + rng.normal(0, 0.8)
                p["x"], p["y"] = x, y
                if R - np.hypot(x, y) < 3 and p["esc"] < -5:                     # entry: into the middle of the nearest gap
                    ph = np.sort([q["phi"] % (2 * np.pi) for q in P_ if q["st"] == "w"]); a0 = np.arctan2(y, x) % (2 * np.pi)
                    if len(ph): ext = np.r_[ph, ph[0] + 2 * np.pi]; i = np.searchsorted(ext, a0) if a0 >= ph[0] else len(ph); i = min(max(i, 1), len(ext) - 1); a0 = 0.5 * (ext[i - 1] + ext[i])
                    p.update(st="w", phi=a0, d=1.5)
            tr.append((p["id"], x, y))
            if rng.random() > 0.02: rows.append((np.round(x * 2) / 2, np.round(y * 2) / 2, 0.5, 0.0, p["id"]))
        truth.append(np.array(tr)); frames.append(np.array(rows) if rows else np.zeros((0, 5)))
    return truth, frames
dist = lambda xy: R - np.hypot(xy[..., 0], xy[..., 1])
def p_esc_truth(truth):
    st = es = 0
    for n in range(0, len(truth) - 10, 5):
        a = {int(i): (x, y) for i, x, y in truth[n]}; b = {int(i): (x, y) for i, x, y in truth[n + 10]}
        for i, xy in a.items():
            if dist(np.array(xy)) < 4: st += 1; es += i in b and dist(np.array(b[i])) > 6
    return es / max(st, 1)
def p_esc_tracks(T, nf):
    st = es = 0
    for tr in T:
        pos = dict(zip(tr["n"], tr["xy"]))
        for n, xy in pos.items():
            if n % 5 or n + 10 >= nf or dist(np.array(xy)) >= 4: continue
            st += 1; e = pos.get(n + 10); es += e is not None and dist(np.array(e)) > 6
    return es / max(st, 1)
def p_esc_chain(fr):
    st = es = 0
    for n in range(0, len(fr) - 10, 5):
        f = fr[n]
        for x, y, q, a, _ in f:
            if dist(np.array([x, y])) < 4:
                st += 1; e = chain_follow([g[:, :4] for g in fr], n, (x, y), 1); es += e is not None and dist(e) > 6
    return es / max(st, 1)
def switches(T, fr):
    lab = [{(round(r_[0], 2), round(r_[1], 2)): int(r_[4]) for r_ in f} for f in fr]; sw = {"bulk": [0, 0], "wall": [0, 0]}; frag = {}
    for tr in T:
        ids = [lab[n].get((round(x, 2), round(y, 2))) for n, (x, y) in zip(tr["n"], tr["xy"])]
        for k in range(1, len(ids)):
            cls = "wall" if dist(np.array(tr["xy"][k])) < 4 else "bulk"; sw[cls][1] += 1; sw[cls][0] += ids[k] != ids[k - 1]
        for i in set(ids): frag[i] = frag.get(i, 0) + 1
    return {k: 1000 * v[0] / max(v[1], 1) for k, v in sw.items()}, np.mean(list(frag.values()))
if __name__ == "__main__":
    L = ["SYNTHETIC GROUND-TRUTH TEST of the defect linkers (synth.py; 3 seeds per scenario, 2500 frames each)",
         "scenario | p_esc truth | p_esc NEW(v2) / NN / CHAIN (ratio to truth) | ID switches per 1000 links, bulk: NEW / NN | wall: NEW / NN | tracks per true defect: NEW / NN"]
    for nm, Nw, r, sf in (("sparse N_w 8", 8, 0.015, "sparse"), ("crowded N_w 28", 28, 0.0012, "crowded")):
        steps = np.load(f"steps_bulk_{sf}.npy"); acc = []
        for seed in (1, 2, 3):
            truth, fr = simulate(Nw, r, seed, steps=steps); frs = [f[:, :4] if len(f) else np.zeros((0, 4)) for f in fr]
            Tn = track_new(frs); To = track_nn(frs); T1 = track_new(frs, dual=False); pt = p_esc_truth(truth)
            pn, po, pc = p_esc_tracks(Tn, len(fr)), p_esc_tracks(To, len(fr)), p_esc_chain(fr)
            (swn, fn), (swo, fo), (sw1, f1) = switches(Tn, fr), switches(To, fr), switches(T1, fr); p1 = p_esc_tracks(T1, len(fr))
            acc.append((pt, pn, po, pc, swn["bulk"], swo["bulk"], swn["wall"], swo["wall"], fn, fo, p1, sw1["bulk"], sw1["wall"]))
        a = np.array(acc).mean(0)
        L.append(f"{nm:15s} | {a[0]:.4f} | {a[1]:.4f} / {a[2]:.4f} / {a[3]:.4f} ({a[1] / a[0]:.2f} / {a[2] / a[0]:.2f} / {a[3] / a[0]:.2f}) | {a[4]:.1f} / {a[5]:.1f} | {a[6]:.1f} / {a[7]:.1f} | {a[8]:.1f} / {a[9]:.1f}"
                 f"  || NEW v1 (prediction only): p_esc ratio {a[10] / a[0]:.2f}, switches bulk {a[11]:.1f}, wall {a[12]:.1f}")
    open("synth.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
