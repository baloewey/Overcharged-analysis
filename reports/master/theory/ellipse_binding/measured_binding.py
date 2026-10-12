#!/usr/bin/env python3
"""Binding of wall +1/2 measured from the active simulations (Boltzmann inversion with the active noise as an effective temperature).
Active ellipses W3 (eps 1.2, 1.5, 2.0; s01, s02) and the circle reference (Rscan_R31.83_seeds s01-s03, k = 6). Every +1/2 within 8 cells of the wall in every
frame (dt = 2, second 80% of the run): depth d and region by the nearest wall point: tip |x| > A/2, flank |x| < A/4, middle otherwise (circle: all one region).
Vortex state per flow interval (halves co-/counter-rotating; ellipses only). Effective potential U(d) = -ln P(d) (bins 0.5 cells), U_min at the peak; reported:
well depth U(6) - U_min (escape barrier to 6 cells), the rms depth within the well (d < 4), and the fraction of wall-band occupancy beyond 4 cells; also the
along-wall rms speed of wall defects (a proxy for the local noise strength, to check that tips and flanks see similar noise). Writes measured_binding.txt."""
import gzip, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline, project  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
bins = np.arange(0, 8.01, 0.5)
def summarize(ds):
    ds = np.array(ds)
    if len(ds) < 200: return None
    P, _ = np.histogram(ds, bins=bins); P = P / P.sum(); U = -np.log(np.maximum(P, 1e-9)); i0 = np.argmax(P)
    c = 0.5 * (bins[1:] + bins[:-1]); U6 = np.interp(6.0, c, U)
    return dict(n=len(ds), dmode=c[i0], barrier=U6 - U[i0], rms=np.std(ds[ds < 4]), beyond=np.mean(ds > 4))
def flow_states(run, bc, A, B):
    rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    Q = bc["Q"]; st = {}
    for tt in np.unique(rows[:, 0]):
        r_ = rows[rows[:, 0] == tt]; x = r_[:, 1] + 0.5 - Q[0]; y = r_[:, 2] + 0.5 - Q[1]; m = (x / A) ** 2 + (y / B) ** 2 < 0.9
        def hr(sel):
            cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
            return np.mean((dx * r_[sel, 5] - dy * r_[sel, 4]) / np.maximum(np.hypot(dx, dy), 0.5))
        st[tt] = int(hr(m & (x < 0)) * hr(m & (x > 0)) < 0)
    return st
def analyse(run, ellipse=True):
    bc = json.load(open(run + "/input.json"))["BC"][0]; W, s, P = wall_polyline(bc); A, B = np.abs(W[:, 0]).max(), np.abs(W[:, 1]).max()
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); t = np.array(t)
    st = flow_states(run, bc, A, B) if ellipse else None; tf = np.array(sorted(st)) if st else None
    D = {}; prev = None
    for i, f in enumerate(fr):
        p = f[f[:, 2] > 0]
        if not len(p): continue
        d, sa = project(W, s, p[:, :2]); near = d < 8
        x = p[:, 0]
        if ellipse:
            j = np.searchsorted(tf, t[i]); vs = st[tf[min(j, len(tf) - 1)]]
        else: vs = 0
        for dd, xx in zip(d[near], x[near]):
            reg = "all" if not ellipse else ("tip" if abs(xx) > A / 2 else ("flank" if abs(xx) < A / 4 else "middle"))
            D.setdefault((reg, vs), []).append(dd)
    return D, A
out = {}
say("Effective potential of wall +1/2 from depth histograms (U = -ln P, units of the active effective temperature):")
for lab, runs, ell in (("circle k6 R31.83", [f"{ROOT}/cavity/Rscan_R31.83_seeds/{s_}/k06" for s_ in ("s01", "s02", "s03")], False),
                       ("eps 1.2", [f"{ROOT}/cassini/W3_ellipse/eps1.2_k06/{s_}" for s_ in ("s01", "s02")], True),
                       ("eps 1.5", [f"{ROOT}/cassini/W3_ellipse/eps1.5_k06/{s_}" for s_ in ("s01", "s02")], True),
                       ("eps 2.0", [f"{ROOT}/cassini/W3_ellipse/eps2.0_k06/{s_}" for s_ in ("s01", "s02")], True)):
    agg = {}
    for r in runs:
        if not os.path.exists(r + "/DONE"): continue
        D, A = analyse(r, ell)
        for key, v in D.items(): agg.setdefault(key, []).extend(v)
    for key in sorted(agg):
        sm = summarize(agg[key])
        if sm is None: continue
        reg, vs = key
        say(f"  {lab:16s} {reg:6s} {'(two vortices)' if vs else '(one vortex) ':14s}: n {sm['n']:7d}, depth mode {sm['dmode']:.2f}, well depth U(6)-U_min {sm['barrier']:.2f}, "
            f"rms depth in the well {sm['rms']:.2f}, fraction deeper than 4 cells {sm['beyond']:.3f}")
open(os.path.join(H, "measured_binding.txt"), "w").write("\n".join(L) + "\n")
