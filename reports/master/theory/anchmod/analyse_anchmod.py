#!/usr/bin/env python3
"""Modulated-winding test (branch anchmod, anchoring angle k phi + A sin(m phi) on the circle, eMode 2).
Prediction (active): the wall +1/2 sit where the anchoring director makes the angle psi_rel(phi) = (k-1) phi + A sin(m phi) with the normal
equal to tau/2 mod pi (tails at the wall shifted by the tilt lock; tau = measured mean tilt, sign = rotation sense); for A = 0 these are
the evenly spaced homeotropic points. Measured: frames with 2k-2 = 10 wall +1/2 (passive: 2k = 12), second half of the run; sites =
circular k-means of the wall positions (10 or 12 clusters, initialised at the predicted or evenly spaced angles); per run: site angles,
rms distance to the predicted sites (cells along the wall), and the spacing pattern. Across initial conditions: rms difference of the site sets.
Writes anchmod_results.txt and anchmod_sites.png."""
import glob, json, os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
def predicted(k, A, m, tau, n_sites):
    f = lambda ph, n: (k - 1) * ph + A * np.sin(m * ph) - n * np.pi - tau / 2
    out = []
    for n in range(-1, n_sites + 2):
        g = lambda ph: f(ph, n)
        grid = np.linspace(-0.5, 2 * np.pi + 0.5, 4000); v = g(grid); idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
        for i in idx: out.append(brentq(g, grid[i], grid[i + 1]) % (2 * np.pi))
    out = np.unique(np.round(np.sort(out), 6)); return out[:n_sites] if len(out) >= n_sites else out
def sites_from(phis, init):
    c = np.array(init, float)
    for _ in range(50):
        d = np.angle(np.exp(1j * (phis[:, None] - c[None, :]))); lab = np.argmin(np.abs(d), 1)
        c = np.array([np.angle(np.mean(np.exp(1j * phis[lab == j]))) % (2 * np.pi) if (lab == j).any() else c[j] for j in range(len(c))])
    return np.sort(c)
def match_rms(a, b, R):
    a, b = np.sort(a), np.sort(b); best = np.inf
    for s in range(len(b)):
        d = np.angle(np.exp(1j * (a - np.roll(b, s)))); best = min(best, np.sqrt(np.mean(d ** 2)))
    return best * R
res = {}
fig, ax = plt.subplots(1, 1, figsize=(9, 4.5)); yrow = 0
for run in sorted(glob.glob(ROOT + "/anchmod_test/*")):
    name = os.path.basename(run)
    if not os.path.exists(run + "/DONE"): say(f"{name}: not finished"); continue
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = int(bc["customAnch"]); A = bc["anchAmp"]; m = int(bc["anchMode"])
    active = name.startswith("active"); n_sites = 2 * k - 2 if active else 2 * k
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.5); P, TL = [], []
    for f in fr:
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); w = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if w.sum() != n_sites: continue
        P.append(ph[w] % (2 * np.pi)); TL.append(np.angle(np.mean(np.exp(1j * (p[w, 3] - ph[w])))))
    if len(P) < 20: say(f"{name}: only {len(P)} frames with {n_sites} wall defects"); continue
    phis = np.concatenate(P); tau = np.angle(np.mean(np.exp(1j * np.array(TL))))
    pred = predicted(k, A, m, tau if active else 0.0, n_sites)
    init = pred if len(pred) == n_sites else (np.arange(n_sites) + 0.5) * 2 * np.pi / n_sites
    s = sites_from(phis, init); spread = np.sqrt(np.mean(np.angle(np.exp(1j * (phis[:, None] - s[None, :])))[np.arange(len(phis)), np.argmin(np.abs(np.angle(np.exp(1j * (phis[:, None] - s[None, :])))), 1)] ** 2)) * R
    rms_pred = match_rms(s, pred, R) if len(pred) == n_sites else np.nan
    rms_even = min(match_rms(s, (np.arange(n_sites) + off) * 2 * np.pi / n_sites, R) for off in np.linspace(0, 1, 50))
    gaps = np.diff(np.append(s, s[0] + 2 * np.pi)) * R
    res[name] = dict(s=s, active=active, A=A, m=m, R=R)
    say(f"{name:28s}: {len(P)} frames, tilt {np.degrees(tau):+.0f} deg; site scatter {spread:.2f} cells; rms to predicted sites {rms_pred:.2f} cells; "
        f"to the best evenly spaced lattice {rms_even:.2f} cells; gaps min/max {gaps.min():.1f}/{gaps.max():.1f} cells")
    ax.plot(np.degrees(s), np.full(len(s), yrow), "o", color="C0" if active else "C2", ms=5)
    if len(pred) == n_sites: ax.plot(np.degrees(pred), np.full(len(pred), yrow), "|", color="C3", ms=12)
    ax.text(362, yrow, name, fontsize=7, va="center"); yrow += 1
ax.set_xlim(0, 360); ax.set_xlabel("site angle (deg)"); ax.set_yticks([]); ax.set_title("measured sites (dots) and predicted homeotropic points shifted by the tilt (red ticks)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(H, "anchmod_sites.png"), dpi=150)
for grp in ("active_A0.5_m1", "active_A1.0_m1", "active_A1.0_m2", "passive_A1.0_m1"):
    names = [n for n in res if n.startswith(grp)]
    if len(names) > 1:
        d = [match_rms(res[a]["s"], res[b]["s"], res[a]["R"]) for i, a in enumerate(names) for b in names[i + 1:]]
        say(f"{grp}: site sets across initial conditions differ by rms {np.mean(d):.2f} cells (pairs: " + ", ".join(f"{x:.2f}" for x in d) + ")")
open(os.path.join(H, "anchmod_results.txt"), "w").write("\n".join(L) + "\n")
