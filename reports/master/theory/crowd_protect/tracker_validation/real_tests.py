#!/usr/bin/env python3
"""Validation of the linkers on real data (2026-10-22). Representative runs: circles and aspect-2 ellipses at rho 0.04 / 0.09 / 0.14 (shape sweep, s01),
channels W24 mixed and sym (channel_current s3).
(1) EXACT CLOSURE. The charge enclosed by a closed curve changes only by crossings (pair events are neutral; a pair born or annihilated straddling the curve
    counts). Curves: d = 4 and d = 8 from the wall (cavities), the mid-line (channels). Tracked estimate: cumulative net charge carried across the curve by
    track crossings + straddling pairs (births/deaths of a +1/2 and a -1/2 within 5 cells, same frame, on opposite sides). Score: rms and final difference
    between tracked and exact over the run (charge units).
(2) FRAME DECIMATION. p_esc (shape_runs definition) from tracks on all frames vs on every 2nd frame (4 t.u.; same 20 t.u. horizon). Robust = ratio ~1.
(3) GATE. p_esc with NEW gates 6, 9, 12.
Linkers: NEW v2 (gate 9), NEW v1 (prediction only), NN (3.5). -> real_tests.txt"""
import json, sys, numpy as np
from multiprocessing import Pool
from trackers import track_new, track_nn
sys.path.insert(0, "../exchange"); import exchange_count as X  # noqa
RUNS = [("circle rho0.04", "shapes_server/circle/circle/Req1.5/k04/s01"), ("circle rho0.14", "shapes_server/circle/circle/Req1.5/k13/s01"),
        ("ellipse2 rho0.04", "shapes_server/ellipse/eps2.00/Req1.5/k05/s01"), ("ellipse2 rho0.09", "shapes_server/ellipse/eps2.00/Req1.5/k10/s01"),
        ("ellipse2 rho0.14", "shapes_server/ellipse/eps2.00/Req1.5/k16/s01"), ("channel W24 mixed", "channel_current/W24_mixed_s3"), ("channel W24 sym", "channel_current/W24_sym_s3")]
def closure(T, fr, side, nf):
    cum = np.zeros(nf)
    for tr in T:
        sd = side(np.array(tr["xy"]))
        for i in np.flatnonzero(sd[1:] != sd[:-1]): cum[tr["n"][i + 1]] += 0.5 * tr["sg"] * (1 if sd[i + 1] else -1)
    for idx in (0, -1):
        ev = {}
        for tr in T:
            n = tr["n"][idx]
            if (idx == 0 and n == 0) or (idx == -1 and n >= nf - 3): continue
            ev.setdefault(n, []).append((tr["sg"], np.array(tr["xy"][idx])))
        for n, es in ev.items():
            Pp = [e for e in es if e[0] > 0]; M = [e for e in es if e[0] < 0]; used = set()
            for _, a in Pp:
                best = None
                for j, (_, b) in enumerate(M):
                    if j not in used and np.hypot(*(a - b)) < 5 and (best is None or np.hypot(*(a - b)) < best[0]): best = (np.hypot(*(a - b)), j, b)
                if best:
                    used.add(best[1]); sa, sb = side(a[None])[0], side(best[2][None])[0]
                    if sa != sb: cum[n if idx == 0 else min(n + 1, nf - 1)] += (0.5 if sa else -0.5) * (1 if idx == 0 else -1)
    return np.cumsum(cum)
def pesc(T, dist, nf, step=1):
    st = es = 0; h = 10 // step; s5 = max(5 // step, 1)
    for tr in T:
        pos = dict(zip(tr["n"], tr["xy"]))
        for n, xy in pos.items():
            if n % s5 or n + h >= nf or dist(np.array([xy]))[0] >= 4: continue
            st += 1; e = pos.get(n + h); es += e is not None and dist(np.array([e]))[0] > 6
    return es / max(st, 1)
def one(item):
    nm, rel = item; run = X.ROOT + "/" + rel; kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nf = len(fr)
    L = json.load(open(run + "/input.json"))["domain"][0] if kind == "channel" else None
    if kind == "channel":
        J = json.load(open(run + "/input.json")); ym = 0.5 * (J["BC"][2]["Q"][1] + J["BC"][3]["Q"][1]); curves = {"mid-line": lambda xy: xy[:, 1] > ym}
    else: curves = {"d=4": lambda xy: dist(xy) > 4, "d=8": lambda xy: dist(xy) > 8}
    lk = {"NEW v2": lambda f: track_new(f, 9.0, L=L), "NEW v1": lambda f: track_new(f, 9.0, L=L, dual=False), "NN": lambda f: track_nn(f, 3.5, L=L)}
    out = {}
    for ln, fn in lk.items():
        T = fn(fr); res = {}
        for cn, side in curves.items():
            ex = np.array([0.5 * np.sign(f[side(f[:, :2]), 2]).sum() if len(f) else np.nan for f in fr]); ex = ex - np.nanmean(ex[:10]); tc = closure(T, fr, side, nf)
            res[cn] = (np.sqrt(np.nanmean((tc - ex) ** 2)), tc[-1] - np.nanmean(ex[-10:]))
        res["p_esc"] = pesc(T, dist, nf) if kind == "cavity" else np.nan
        if kind == "cavity":
            T2 = fn(fr[::2]); T2 = [dict(tr, n=[2 * n for n in tr["n"]]) for tr in T2]; res["p_esc_4tu"] = pesc(T2, dist, nf, step=1)
        out[ln] = res
    if kind == "cavity":
        out["gates"] = {g: pesc(track_new(fr, g), dist, nf) for g in (6.0, 9.0, 12.0)}
    return nm, out
if __name__ == "__main__":
    with Pool(7) as p: R = p.map(one, RUNS)
    Lt = ["REAL-DATA VALIDATION (real_tests.py)", "", "(1) EXACT CLOSURE: tracked minus exact net charge carried across the curve: rms over the run / final (charge units; a perfect linker gives ~0)"]
    for nm, o in R:
        for cn in [k for k in o["NN"] if k.startswith("d=") or k == "mid-line"]:
            Lt.append(f"  {nm:18s} {cn:8s}: " + " | ".join(f"{ln} {o[ln][cn][0]:5.1f} / {o[ln][cn][1]:+6.1f}" for ln in ("NEW v2", "NEW v1", "NN")))
    Lt.append("\n(2) FRAME DECIMATION: p_esc from all frames (2 t.u.) -> every 2nd frame (4 t.u.); ratio")
    for nm, o in R:
        if "p_esc_4tu" in o["NN"]: Lt.append(f"  {nm:18s}: " + " | ".join(f"{ln} {o[ln]['p_esc']:.4f} -> {o[ln]['p_esc_4tu']:.4f} (x{o[ln]['p_esc_4tu'] / max(o[ln]['p_esc'], 1e-9):.2f})" for ln in ("NEW v2", "NEW v1", "NN")))
    Lt.append("\n(3) GATE (NEW v2): p_esc at gates 6 / 9 / 12 cells")
    for nm, o in R:
        if "gates" in o: Lt.append(f"  {nm:18s}: " + " / ".join(f"{v:.4f}" for v in o["gates"].values()))
    open("real_tests.txt", "w").write("\n".join(Lt) + "\n"); print("\n".join(Lt))
