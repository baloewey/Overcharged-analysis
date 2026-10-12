#!/usr/bin/env python3
"""Stage 2, combined analysis of all feedback-matrix runs (first round M_e*.json, enlarged ensembles M_2_c*.json, realistic disorder
M_3_c*.json, measured snapshots s01 and s02). Same mode classification as analyse_stage2.py. Per eps and disorder class:
P(top mode two-vortex) with a 95% Wilson interval, gain ratio one/two (median, 10-90%). Measured snapshots: per seed, as in analyse_stage2.py.
Also a figure stage2.png: P(two-vortex top) and median gain ratio vs sigma_t per eps. Writes stage2b_summary.txt."""
import glob, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
ns = {"__file__": os.path.join(HERE, "analyse_stage2.py")}
exec(open(os.path.join(HERE, "analyse_stage2.py")).read().split("for tag in (")[0], ns)
modes, halves, read_with_angle, RUNS = ns["modes"], ns["halves"], ns["read_with_angle"], ns["RUNS"]
L = []
def say(s): print(s); L.append(s)
def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0, c - h), min(1, c + h))
EPS = {"e10": 1.0, "e12": 1.2, "e15": 1.5, "e20": 2.0, "e25": 2.5}
recs = []                                                    # (tag, kind, sigma, vac, seed, interval, t, mode-dict, sites)
for f in sorted(glob.glob(os.path.join(HERE, "M_*.json"))):
    name = os.path.basename(f)[2:-5]
    cl = json.load(open(os.path.join(HERE, f"cfg_{name}.json")))    # ids repeat across eps inside a chunk: match by order
    for c, r in zip(cl, json.load(open(f))):
        assert c["id"] == r["id"]; tag = c.get("tag", name)
        x = np.array([s[0] for s in c["sites"]])
        recs.append(dict(tag=tag, kind=c["kind"], sigma=c.get("sigma", 1.0 if c["kind"] == "vacancy" else 0.0), vac=c.get("vac", 0),
                         seed=c.get("seed", "s01"), interval=c.get("interval"), t=c.get("t"), m=modes(np.array(r["M"]), x), sites=c["sites"]))
say(f"{len(recs)} configurations")
groups = {}
for r in recs:
    if r["kind"] in ("disorder", "vacancy", "realistic", "perfect"):
        groups.setdefault((r["tag"], r["sigma"], r["vac"]), []).append(r["m"])
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for tag in EPS:
    keys = sorted(k for k in groups if k[0] == tag)
    if not keys: continue
    say(f"\n=== eps {EPS[tag]}")
    xs, ps, rs = [], [], []
    for key in keys:
        ms = groups[key]; n = len(ms); k2 = sum(m["typ"][0] == "two" for m in ms); rat = np.array([m["gu"] / m["gd"] for m in ms])
        lo, hi = wilson(k2, n)
        say(f"  sigma_t {key[1]:>3g}, vacancies {key[2]}: n {n:3d}, P(top two-vortex) {k2 / n:.2f} [{lo:.2f}-{hi:.2f}], "
            f"gain ratio one/two median {np.median(rat):.3f} (10-90% {np.percentile(rat, 10):.3f}-{np.percentile(rat, 90):.3f}), frac ratio<1 {np.mean(rat < 1):.2f}")
        if key[2] == 0: xs.append(key[1]); ps.append(k2 / n); rs.append(np.median(rat))
    ax[0].plot(xs, ps, "o-", label=f"$\\varepsilon$ = {EPS[tag]}"); ax[1].plot(xs, rs, "o-", label=f"$\\varepsilon$ = {EPS[tag]}")
ax[0].set_xlabel(r"tangential disorder $\sigma_t$ (cells)"); ax[0].set_ylabel("P(top mode two-vortex)"); ax[0].set_title("(a) imperfect lattices, no vacancies", fontsize=10)
ax[1].set_xlabel(r"$\sigma_t$ (cells)"); ax[1].set_ylabel("gain ratio one/two (median)"); ax[1].axhline(1, color="k", lw=0.6); ax[1].set_title("(b) one-vortex margin", fontsize=10)
ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(HERE, "stage2.png"), dpi=150)
SF = os.path.join(HERE, "..", "..", "stokes_frames")
for seed in ("s01", "s02"):
    for tag in ("e12", "e15", "e20"):
        ms = sorted([r for r in recs if r["kind"] == "measured" and r["seed"] == seed and r["tag"] == tag], key=lambda r: r["interval"])
        if not ms: continue
        if seed == "s01": Hh = halves(tag)
        else:
            d = np.load(os.path.join(HERE, f"s02_{tag}_vel.npz")); x, y = d["x"], d["y"]; Hh = []
            for vx, vy in zip(d["vx"], d["vy"]):
                rr = []
                for msk in (x < 0, x > 0):
                    cx, cy = x[msk].mean(), y[msk].mean(); dx, dy = x[msk] - cx, y[msk] - cy
                    rr.append(np.mean((dx * vy[msk] - dy * vx[msk]) / np.maximum(np.hypot(dx, dy), 0.5)))
                Hh.append(rr)
            Hh = np.array(Hh)
        co = np.array([Hh[r["interval"], 0] * Hh[r["interval"], 1] > 0 for r in ms]); two = np.array([r["m"]["typ"][0] == "two" for r in ms])
        rat = np.array([r["m"]["gu"] / r["m"]["gd"] for r in ms])
        say(f"\nmeasured {seed} eps {EPS[tag]}: {len(ms)} snapshots, wall +1/2 {np.mean([len(r['sites']) for r in ms]):.1f}; measured one-vortex in {co.mean():.0%} "
            f"of intervals; predicted top two-vortex in {two.mean():.0%}; gain ratio median {np.median(rat):.3f} "
            f"(one-vortex intervals {np.median(rat[co]) if co.any() else np.nan:.3f}, two-vortex {np.median(rat[~co]) if (~co).any() else np.nan:.3f})")
open(os.path.join(HERE, "stage2b_summary.txt"), "w").write("\n".join(L) + "\n")
