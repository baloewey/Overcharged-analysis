#!/usr/bin/env python3
"""
Wall-layer and docking analysis for circular cavities with imposed winding k (H6, H7).

Per active run directory (input.json + defects.dat[.gz]):

  Edge/bulk division, distance from the wall s = R - r:
    passive  d_p  = 95th percentile of s over defects in the passive run of the same (R, k)
                    (--passive DIR); "no layer" if fewer than half of the passive defects lie
                    within --layer-max cells of the wall
    density  d_ρ  = minimum of the radial density of same-sign defects (sign of k) between the
                    wall peak and the bulk (data-driven, needs no passive run)
  Bulk charge q_b(t) (net charge of defects with s > d), with block-average errors (10 blocks):
    <q_b>, P(q_b = +1)
  Layer (H6): occupancy N_L(t) (same-sign defects with s <= d), mean, relative fluctuation
    std/mean, radial width (std of s inside the layer); compared with the passive run
  Trajectories: frame-to-frame linking of same-sign defects (optimal assignment, max jump
    --rlink cells); unmatched = birth/death (pair creation, annihilation, tracker misses)
  Docking (H7), for defects of the sign of k:
    dock = a track crossing from the bulk into the layer, undock = the reverse; rates per unit time
    circulation time = time a docking track spent in the edge zone d < s <= d + --zone before docking
    azimuthal speed |dθ/dt| r of bulk defects vs s (binned)
    gap test: at each dock, rank of the gap (between neighbouring layer defects in angle) that
    the docking defect enters; fraction docking into the largest gap vs random expectation 1/n_gaps

Usage:
  layer_analysis.py RUN [RUN ...] [--passive-root DIR] [--tfrac 0.2] [--out layer.csv]
    --passive-root: directory holding passive runs as <root>/kNN (or <root>/rhoX/kNN/s01 for sweep
    layouts: the same relative rho/k path is used). Without it only the density edge is used.
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from defect_audit import read_frames  # noqa: E402


def read_input(rundir):
    with open(os.path.join(rundir, "input.json")) as f:
        inp = json.load(f)
    bc = inp["BC"][0]
    ca = bc.get("customAnch", 0) or 1
    return {"C": bc["Q"][:2], "R": bc["R"], "k": ca, "act": inp["species"][0]["act"],
            "dt_out": inp.get("defectsOut", 50) * inp.get("dt", 0.1)}


def load(rundir, tfrac):
    p = read_input(rundir)
    t, fr = read_frames(os.path.join(rundir, "defects.dat"))
    keep = t >= t.min() + tfrac * (t.max() - t.min()) if len(t) else []
    t = t[keep]
    fr = [f for f, kp in zip(fr, keep) if kp]
    C, R = p["C"], p["R"]
    out = []
    for f in fr:
        if len(f):
            x, y = f[:, 0] - C[0], f[:, 1] - C[1]
            out.append(np.c_[x, y, np.hypot(x, y), np.arctan2(y, x), np.sign(f[:, 2])])
        else:
            out.append(np.zeros((0, 5)))
    return p, t, out   # columns: x, y, r, theta, sign


# ------------------------------------------------------------------ edges
def passive_edge(prun, tfrac, layer_max=5.0):
    p, t, fr = load(prun, tfrac)
    sk = np.sign(p["k"])
    s = np.concatenate([p["R"] - f[f[:, 4] == sk, 2] for f in fr if len(f)]) if fr else np.array([])
    if len(s) == 0:
        return np.nan, False
    exists = np.mean(s <= layer_max) >= 0.5
    return float(np.percentile(s, 95)), bool(exists)


def density_edge(p, fr, dbin=0.5):
    R, sk = p["R"], np.sign(p["k"])
    s = np.concatenate([R - f[f[:, 4] == sk, 2] for f in fr if len(f)])
    edges = np.arange(0, min(R, 12) + dbin, dbin)
    h, _ = np.histogram(s, edges)
    rin, rout = R - edges[1:], R - edges[:-1]
    dens = h / (np.pi * (rout**2 - rin**2)) / len(fr)
    mid = 0.5 * (edges[1:] + edges[:-1])
    ipk = int(np.argmax(dens[mid < 3])) if np.any(mid < 3) else 0
    sel = np.where((mid > mid[ipk]) & (mid <= min(8, R / 2)))[0]
    if len(sel) == 0:
        return np.nan
    j = sel[np.argmin(np.convolve(dens, np.ones(3) / 3, "same")[sel])]
    return float(mid[j])


# ------------------------------------------------------------------ stats
def block(x, nb=10):
    x = np.asarray(x, float)
    m = len(x) // nb
    if m < 1:
        return x.mean() if len(x) else np.nan, np.nan
    b = x[: m * nb].reshape(nb, m).mean(1)
    return float(x.mean()), float(b.std(ddof=1) / np.sqrt(nb))


# ------------------------------------------------------------------ linking
def link(t, fr, rlink, max_gap=3):
    """Frame-to-frame linking of same-sign defects by optimal assignment (Hungarian), with gap
    closing: a track may miss up to max_gap frames (tracker flicker at the wall) and still be
    continued from its last position. Returns tracks: arrays (frame index, x, y, r, theta, sign)."""
    tracks = []
    live = {}                     # track id -> (last frame index, last row)
    for n, f in enumerate(fr):
        ids = [i for i, (m, _) in live.items() if n - m <= max_gap + 1]
        used = set()
        if ids and len(f):
            P = np.array([live[i][1] for i in ids])
            d = np.hypot(P[:, None, 0] - f[None, :, 0], P[:, None, 1] - f[None, :, 1])
            gapf = np.array([n - live[i][0] for i in ids])[:, None]
            d[(P[:, None, 4] != f[None, :, 4]) | (d > rlink * np.sqrt(gapf))] = 1e9
            cost = d + 0.5 * (gapf - 1)          # prefer continuing tracks seen in the last frame
            ri, ci = linear_sum_assignment(cost)
            for a_, b_ in zip(ri, ci):
                if d[a_, b_] < 1e9:
                    tracks[ids[a_]].append((n, *f[b_]))
                    live[ids[a_]] = (n, f[b_])
                    used.add(b_)
        for b_ in range(len(f)):
            if b_ not in used:
                tracks.append([(n, *f[b_])])
                live[len(tracks) - 1] = (n, f[b_])
        live = {i: v for i, v in live.items() if n - v[0] <= max_gap}
    return [np.array(tr) for tr in tracks]


def docking(p, t, fr, tracks, d, zone):
    R, sk = p["R"], np.sign(p["k"])
    dt = np.median(np.diff(t)) if len(t) > 1 else np.nan
    docks, undocks, circ, gap_top, gap_n = 0, 0, [], 0, []
    for tr in tracks:
        if len(tr) < 2 or tr[0, 5] != sk:
            continue
        s = R - tr[:, 3]
        inlayer = s <= d
        for i in range(1, len(tr)):
            if inlayer[i] and not inlayer[i - 1]:
                docks += 1
                # time in the edge zone just before docking
                j = i - 1
                while j >= 0 and d < s[j] <= d + zone:
                    j -= 1
                circ.append((i - 1 - j) * dt)
                # gap test: layer defects (same sign) in frame tr[i,0], excluding the docker
                n = int(tr[i, 0])
                f = fr[n]
                L = f[(R - f[:, 2] <= d) & (f[:, 4] == sk)]
                th = np.sort(L[:, 3])
                th_d = tr[i, 4]
                others = th[np.abs(np.angle(np.exp(1j * (th - th_d)))) > 1e-9]
                if len(others) >= 2:
                    ths = np.sort(others)
                    gaps = np.diff(np.r_[ths, ths[0] + 2 * np.pi])
                    # gap containing the docking angle
                    rel = (th_d - ths[0]) % (2 * np.pi)
                    idx = np.searchsorted(np.cumsum(gaps), rel)
                    idx = min(idx, len(gaps) - 1)
                    gap_top += int(gaps[idx] >= gaps.max() - 1e-12)
                    gap_n.append(len(gaps))
            elif inlayer[i - 1] and not inlayer[i]:
                undocks += 1
    T = (t.max() - t.min()) if len(t) > 1 else np.nan
    # azimuthal speed of bulk same-sign defects vs s
    sbins = np.arange(0, min(R, 16) + 2, 2.0)
    vs = [[] for _ in sbins[:-1]]
    for tr in tracks:
        if len(tr) < 2 or tr[0, 5] != sk:
            continue
        s = R - tr[:, 3]
        dth = np.angle(np.exp(1j * np.diff(tr[:, 4])))
        steps = np.diff(tr[:, 0])
        v = np.abs(dth) * 0.5 * (tr[1:, 3] + tr[:-1, 3]) / (steps * dt)
        sm = 0.5 * (s[1:] + s[:-1])
        for vi, si, st in zip(v, sm, steps):
            if st <= 2:
                b = np.searchsorted(sbins, si) - 1
                if 0 <= b < len(vs):
                    vs[b].append(vi)
    vprof = [(float(0.5 * (sbins[b] + sbins[b + 1])), float(np.mean(v)) if v else np.nan, len(v))
             for b, v in enumerate(vs)]
    nd = len(gap_n)
    return {"dock_rate": docks / T if T else np.nan, "undock_rate": undocks / T if T else np.nan,
            "n_docks": docks, "circ_time_mean": float(np.mean(circ)) if circ else np.nan,
            "circ_time_median": float(np.median(circ)) if circ else np.nan,
            "dock_largest_gap_frac": gap_top / nd if nd else np.nan,
            "dock_largest_gap_random": float(np.mean([1 / g for g in gap_n])) if nd else np.nan,
            "vphi_profile": vprof}


# ------------------------------------------------------------------ per run
def analyse(run, passive_dir, tfrac, rlink, zone, layer_max):
    p, t, fr = load(run, tfrac)
    R, sk = p["R"], np.sign(p["k"])
    res = {"run": run.rstrip("/"), "k": p["k"], "R": R, "rho": abs(p["k"]) / (np.pi * R),
           "act": p["act"], "frames": len(fr)}
    res["d_density"] = density_edge(p, fr)
    if passive_dir and os.path.isdir(passive_dir):
        res["d_passive"], res["passive_layer"] = passive_edge(passive_dir, tfrac, layer_max)
        pp, pt, pfr = load(passive_dir, tfrac)
        dP = res["d_passive"]
        nlp = [np.sum((pp["R"] - f[:, 2] <= dP) & (f[:, 4] == sk)) for f in pfr]
        wp = np.concatenate([pp["R"] - f[(pp["R"] - f[:, 2] <= dP) & (f[:, 4] == sk), 2] for f in pfr if len(f)])
        res["passive_NL_mean"] = float(np.mean(nlp))
        res["passive_width"] = float(np.std(wp)) if len(wp) else np.nan
    else:
        res["d_passive"], res["passive_layer"] = np.nan, None
    for tag in ("passive", "density"):
        d = res[f"d_{tag}"]
        if not np.isfinite(d):
            continue
        qb = np.array([0.5 * f[R - f[:, 2] > d, 4].sum() if len(f) else 0.0 for f in fr])
        m, e = block(qb)
        P1, eP1 = block(np.abs(qb - 1) < 0.25)
        NL = np.array([np.sum((R - f[:, 2] <= d) & (f[:, 4] == sk)) for f in fr])
        w = np.concatenate([R - f[(R - f[:, 2] <= d) & (f[:, 4] == sk), 2] for f in fr if len(f)])
        res.update({f"qb_{tag}": m, f"qb_{tag}_err": e, f"P1_{tag}": P1, f"P1_{tag}_err": eP1,
                    f"NL_{tag}": float(NL.mean()), f"NL_relfluct_{tag}": float(NL.std() / NL.mean()) if NL.mean() else np.nan,
                    f"width_{tag}": float(np.std(w)) if len(w) else np.nan})
    d_use = res["d_passive"] if np.isfinite(res.get("d_passive", np.nan)) else res["d_density"]
    tracks = link(t, fr, rlink)
    lengths = np.array([len(tr) for tr in tracks])
    res["dt_out"] = float(np.median(np.diff(t))) if len(t) > 1 else np.nan
    res["linking_ok"] = bool(res["dt_out"] <= 2.5)   # bulk +1/2 move ~0.8 cells per time unit
    res.update({"edge_used": "passive" if d_use == res.get("d_passive") else "density",
                "n_tracks": len(tracks), "track_len_median": float(np.median(lengths)) if len(lengths) else np.nan})
    if np.isfinite(d_use):
        res.update(docking(p, t, fr, tracks, d_use, zone))
    return res


COLS = ["run", "k", "R", "rho", "act", "frames", "d_passive", "passive_layer", "d_density",
        "qb_passive", "qb_passive_err", "P1_passive", "P1_passive_err",
        "qb_density", "qb_density_err", "P1_density", "P1_density_err",
        "NL_passive", "passive_NL_mean", "NL_relfluct_passive", "width_passive", "passive_width",
        "NL_density", "NL_relfluct_density", "width_density", "edge_used", "dt_out", "linking_ok", "n_tracks", "track_len_median",
        "dock_rate", "undock_rate", "n_docks", "circ_time_mean", "circ_time_median",
        "dock_largest_gap_frac", "dock_largest_gap_random"]


def passive_for(run, root):
    """Map an active run dir to its passive counterpart under root (same kNN, or same rho/k path)."""
    if not root:
        return None
    parts = os.path.normpath(run).split(os.sep)
    kpart = [q for q in parts if q.startswith("k") and q[1:].lstrip("m").isdigit()]
    rpart = [q for q in parts if q.startswith("rho")]
    cands = []
    if rpart and kpart:
        cands.append(os.path.join(root, rpart[-1], kpart[-1], "s01"))
    if kpart:
        cands.append(os.path.join(root, kpart[-1]))
    for c in cands:
        if os.path.isdir(c):
            return c
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--passive-root")
    ap.add_argument("--tfrac", type=float, default=0.2, help="drop the first fraction of production frames")
    ap.add_argument("--rlink", type=float, default=3.0, help="max jump between frames for linking (cells)")
    ap.add_argument("--zone", type=float, default=2.0, help="edge-zone width for circulation time (cells)")
    ap.add_argument("--layer-max", type=float, default=5.0)
    ap.add_argument("--out", default="layer_analysis.csv")
    ap.add_argument("--profiles", help="also write azimuthal-speed profiles (JSON) to this file")
    a = ap.parse_args()
    rows, prof = [], {}
    for r in a.runs:
        res = analyse(r, passive_for(r, a.passive_root), a.tfrac, a.rlink, a.zone, a.layer_max)
        prof[res["run"]] = res.pop("vphi_profile", None)
        rows.append(res)
    with open(a.out, "w") as f:
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in COLS) + "\n")
    if a.profiles:
        json.dump(prof, open(a.profiles, "w"), indent=1)
    print(f"{'run':<34} {'k':>3} {'d_p':>5} {'d_ρ':>5} | {'qb(p)':>11} {'P1(p)':>5} | {'qb(ρ)':>11} | "
          f"{'NL':>5} {'fl':>4} {'w':>4} {'wP':>4} | {'dock/t':>7} {'circ':>5} {'gap':>9}")
    for r in rows:
        g = lambda k, fmt="{:.2f}": (fmt.format(r[k]) if isinstance(r.get(k), (int, float)) and np.isfinite(r[k]) else "–")
        print(f"{r['run'][-34:]:<34} {r['k']:>3} {g('d_passive', '{:.1f}'):>5} {g('d_density', '{:.1f}'):>5} | "
              f"{g('qb_passive')}±{g('qb_passive_err')} {g('P1_passive'):>5} | {g('qb_density')}±{g('qb_density_err')} | "
              f"{g('NL_'+r['edge_used'], '{:.1f}'):>5} {g('NL_relfluct_'+r['edge_used']):>4} {g('width_'+r['edge_used'], '{:.1f}'):>4} "
              f"{g('passive_width', '{:.1f}'):>4} | {g('dock_rate', '{:.3f}'):>7} {g('circ_time_median', '{:.0f}'):>5} "
              f"{g('dock_largest_gap_frac')}/{g('dock_largest_gap_random')}")
    if any(not r.get("linking_ok", True) for r in rows):
        print("WARNING: some runs output defects too rarely for reliable linking (dt_out > 2.5): "
              "docking numbers for those are unreliable (use defectsOut <= 20).")
    print(f"\nwritten {a.out}. d in cells from the wall; qb ± block error; NL layer occupancy, fl = std/mean, "
          "w layer width (wP passive); dock/t per MPCD time; circ = median edge-zone time before docking; "
          "gap = fraction docking into the largest gap / random expectation")


if __name__ == "__main__":
    main()
