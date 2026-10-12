#!/usr/bin/env python3
"""Stage 2 analysis: feedback modes of imperfect wall lattices.
For each configuration: Ms = (M + M^T)/2 (Stokes reciprocity; the FEM asymmetry is reported), eigenvalues l1 >= l2 >= ...
Mode type of an eigenvector v (sites i, x_i): overlap with the uniform pattern c_u = |sum v_i| / sqrt(n) and with the left/right
dipolar pattern c_d = |sum sign(x_i) v_i| / sqrt(n) (|v| = 1). Type 'one-vortex' if c_u > c_d, else 'two-vortex'.
Gain margin: g_u = max l over one-vortex modes, g_d = max l over two-vortex modes; ratio g_u/g_d (< 1: two vortices preferred).
Measured snapshots: compared with the measured flow of the same interval (halves co-rotating or counter-rotating, stokes_frames/<tag>_vel.npz)
and with the measured tilt pattern of the wall +1/2 in the snapshot (sin of the head angle from the outward normal, ccw positive).
Writes stage2_summary.txt."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, "..", "..", "stokes_frames")
RUNS = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle  # noqa
RUNREL = {"e12": "cassini/W3_ellipse/eps1.2_k06/s01", "e15": "cassini/W3_ellipse/eps1.5_k06/s01", "e20": "cassini/W3_ellipse/eps2.0_k06/s01"}
L = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); L.append(s)


def modes(M, x):
    asym = np.abs(M - M.T).max() / np.abs(M).max()
    Ms = (M + M.T) / 2; l, V = np.linalg.eigh(Ms); o = np.argsort(l)[::-1]; l, V = l[o], V[:, o]
    n = len(x); sg = np.where(x < 0, 1.0, -1.0)
    cu = np.abs(V.sum(0)) / np.sqrt(n); cd = np.abs(sg @ V) / np.sqrt(n)
    typ = np.where(cu > cd, "one", "two")
    gu = l[typ == "one"].max() if (typ == "one").any() else np.nan
    gd = l[typ == "two"].max() if (typ == "two").any() else np.nan
    return dict(l=l, V=V, typ=typ, cu=cu, cd=cd, gu=gu, gd=gd, asym=asym)


def halves(tag):
    d = np.load(os.path.join(SF, tag + "_vel.npz")); x, y = d["x"], d["y"]; out = []
    for vx, vy in zip(d["vx"], d["vy"]):
        r = []
        for m in (x < 0, x > 0):
            cx, cy = x[m].mean(), y[m].mean(); dx, dy = x[m] - cx, y[m] - cy
            r.append(np.mean((dx * vy[m] - dy * vx[m]) / np.maximum(np.hypot(dx, dy), 0.5)))
        out.append(r)
    return np.array(out)


for tag in ("e10", "e12", "e15", "e20"):
    f = os.path.join(HERE, f"M_{tag}.json")
    if not os.path.exists(f): continue
    cfg = {c["id"]: c for c in json.load(open(os.path.join(HERE, f"cfg_{tag}.json")))}
    res = json.load(open(f))
    say(f"\n=== {tag} (eps {cfg['perfect']['A'] / cfg['perfect']['B']:.2f}), {len(res)} configurations")
    R = {}
    for r in res:
        c = cfg[r["id"]]; x = np.array([s[0] for s in c["sites"]]); R[r["id"]] = (c, modes(np.array(r["M"]), x))
    say(f"max FEM asymmetry |M - M^T|/max|M|: {max(m['asym'] for c, m in R.values()):.3f}")
    c, m = R["perfect"]
    say(f"perfect lattice: top mode {m['typ'][0]}-vortex, gain ratio one/two {m['gu'] / m['gd']:.3f}")
    groups = {}
    for k, (c, m) in R.items():
        if c["kind"] in ("disorder", "vacancy"):
            g = f"{c['kind']} {c.get('sigma', c.get('vac'))}"
            groups.setdefault(g, []).append(m)
    for g, ms in groups.items():
        top2 = np.mean([m["typ"][0] == "two" for m in ms]); rat = np.array([m["gu"] / m["gd"] for m in ms])
        say(f"{g:14s}: P(top mode two-vortex) = {top2:.2f} ({len(ms)} samples); gain ratio one/two median {np.median(rat):.3f}, "
            f"10-90% {np.percentile(rat, 10):.3f}-{np.percentile(rat, 90):.3f}")
    meas = sorted([(c["interval"], c, m) for c, m in R.values() if c["kind"] == "measured"], key=lambda z: z[0])
    if meas:
        H = halves(tag)
        run = os.path.join(RUNS, RUNREL[tag]); bc = json.load(open(run + "/input.json"))["BC"][0]
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
        rows = []
        for n, c, m in meas:
            co = H[n, 0] * H[n, 1] > 0                                          # measured: halves co-rotating = one vortex
            f = fr[int(np.abs(t - c["t"]).argmin())]
            tilt = []
            for sx, sy, tx, ty in c["sites"]:
                k = np.argmin(np.hypot(f[:, 0] - sx, f[:, 1] - sy)); ang = f[k, 3]
                tilt.append(np.cos(ang) * tx + np.sin(ang) * ty)               # tangential (ccw) component of the head direction
            tilt = np.array(tilt); tn = tilt / max(np.linalg.norm(tilt), 1e-9)
            proj_top = abs(tn @ m["V"][:, 0]) ** 2                             # share of the tilt pattern in the top mode
            rows.append((co, m["typ"][0] == "two", m["gu"] / m["gd"], proj_top, len(c["sites"])))
        rows = np.array(rows, dtype=float)
        co, two, rat, pt, nw = rows.T
        say(f"measured snapshots: {len(rows)}; wall +1/2 per snapshot {nw.mean():.1f}; measured flow one-vortex in {co.mean():.0%} of intervals")
        say(f"   predicted top mode two-vortex in {two.mean():.0%} of snapshots; gain ratio one/two median {np.median(rat):.3f}")
        if 0 < co.mean() < 1:
            say(f"   gain ratio one/two: intervals with one vortex {np.median(rat[co == 1]):.3f}, with two vortices {np.median(rat[co == 0]):.3f}")
            agree = np.mean((two == 1) == (co == 0)); say(f"   predicted vortex number agrees with the measured flow in {agree:.0%} of intervals "
                                                          f"(chance {max(co.mean(), 1 - co.mean()):.0%} by always guessing the majority)")
        say(f"   share of the measured tilt pattern in the predicted top mode: median {np.median(pt):.2f} (random pattern 1/n = {np.mean(1 / nw):.2f})")
open(os.path.join(HERE, "stage2_summary.txt"), "w").write("\n".join(L) + "\n")
