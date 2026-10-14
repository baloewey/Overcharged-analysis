#!/usr/bin/env python3
"""Fig. cont_ellipse: active continuum model (Q-tensor, backflow, one constant, k 6, lambda 2) in ellipses vs the MPCD ellipse sweep at R_eq 1.5 l_a (medians).
Recovered 2026-10-13 from the session transcript (the figure was made by an inline script on 2026-10-09)."""
import csv, glob, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
R = {}
for f in glob.glob("/home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt/ellipse_runs/eps*_z*.json"):
    J = json.load(open(f)); a = J["args"]; o = J["obs"][len(J["obs"]) // 2:]
    R[(a["eps"], a["zeta"])] = (10 - np.mean([x["nwall"] for x in o]), np.nanmean([x.get("coh", np.nan) for x in o]), np.mean([x["two"] for x in o]), np.mean([abs(x["circ"]) for x in o]))
rows = list(csv.DictReader(open("/home/baloewe/Physics/Overcharged/reports/master/theory/shapes/shape_runs.csv")))
def med(key, fam="ellipse"):
    out = {}
    for r in rows:
        if r["family"] in (fam, "circle") and float(r["Req_la"]) == 1.5 and r[key] not in ("", "nan"): out.setdefault(float(r["aspect"]), []).append(float(r[key]))
    xs = sorted(out); return xs, [np.median(out[x]) for x in xs]
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.9))
for i, (lab, j, key) in enumerate((("vacancies ($2k-2$ minus wall count)", 0, "vac"), ("tilt coherence", 1, "coh"), ("P(two vortices)", 2, "two"))):
    for z, mk, col in ((0.3, "o-", "#1f5fa8"), (1.0, "s-", "#c2452d")):
        e = sorted(x for x, zz in R if zz == z); ax[i].plot(e, [R[(x, z)][j] for x in e], mk, color=col, ms=4, label=f"continuum, $\\zeta={z}$")
    xs, ys = med(key); ax[i].plot(xs, ys, "^--", color="0.45", ms=4, label="MPCD ellipses, $R_{eq}=1.5\\,\\ell_a$ (median)")
    ax[i].set_xlabel("aspect ratio"); ax[i].set_title(f"({'abc'[i]}) {lab}", fontsize=9)
fig.legend(*ax[0].get_legend_handles_labels(), fontsize=7, frameon=False, ncol=3, loc="lower center"); fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig("/home/baloewe/Physics/Overcharged/reports/theory_report/figs/fig_cont_ellipse.png", dpi=200)
