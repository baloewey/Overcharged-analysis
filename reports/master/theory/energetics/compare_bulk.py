#!/usr/bin/env python3
"""Depleting escapes: activity everywhere (fineflow + fineflow2, tag _ff12) vs bulk-only activity (wallActOff 4, runs/bulkfine, tag _bulk), from selfflow{tag}.txt and
shear_trigger{tag}.txt. Panels: (a) cos tau and depth; (b) inward flow at the defect and in the ring; (c) defect inward velocity and own-jet estimate; (d) local wall flow
(signed by the cavity circulation), escaper vs control. -> fig_compare_bulk.png"""
import os, re, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
def sf(tag):
    L = open(os.path.join(H, f"selfflow{tag}.txt")).read().splitlines(); jet = float(re.search(r"median ([-\d.]+)", L[0]).group(1))
    R = np.array([[float(x) for x in l.split()] for l in L if re.match(r"^\s*[-+]?\d+\s", l)]); return jet, R
def st(tag):
    R = []
    for l in open(os.path.join(H, f"shear_trigger{tag}.txt")).read().splitlines()[1:]:
        v = [float(x) for x in re.findall(r"[-+]?(?:nan|\d*\.\d+|\d+)", l.replace("±", " "))]; R.append(v)
    return np.array(R)
sets = (("_ff12", "activity everywhere", "#c2452d"), ("_bulk", "passive wall band", "#1f5fa8"))
plt.rcParams.update({"font.size": 11}); fig, ax = plt.subplots(2, 2, figsize=(13, 9.5)); ax = ax.ravel()
for tag, nm, c in sets:
    if not os.path.exists(os.path.join(H, f"selfflow{tag}.txt")): continue
    jet, R = sf(tag); t = R[:, 0]
    ax[0].plot(t, R[:, 2], "-o", color=c, ms=3, label=f"{nm.split(' (')[0]}: cos τ"); a0b = ax[0].twinx() if tag == "_ff12" else a0b; a0b.plot(t, R[:, 1], ":", color=c, label="depth")
    ax[1].plot(t, R[:, 3], "-o", color=c, ms=3, label=f"{nm.split(' (')[0]}: at the defect"); ax[1].plot(t, R[:, 4], "--s", color=c, ms=3, label="ring 4–8 cells")
    ax[2].plot(t, R[:, 6], "-^", color=c, ms=3, label=f"{nm.split(' (')[0]}: defect"); ax[2].plot(t, R[:, 7], "--", color=c, label=f"own-jet estimate ({jet:.2f}·cos τ)")
    if os.path.exists(os.path.join(H, f"shear_trigger{tag}.txt")):
        S = st(tag); ax[3].plot(S[:, 0], S[:, 1], "-o", color=c, ms=3, label=f"{nm.split(' (')[0]}: escape site"); ax[3].plot(S[:, 0], S[:, 4], "--", color=c, label="control site")
a0b.set_ylabel("depth (cells, dotted)")
for a, yl, ti in zip(ax, ("cos τ (1 = head into the bulk)", "inward flow (cells/t.u.)", "inward velocity (cells/t.u.)", "wall flow (cells/t.u.)"),
                     ("(a) tilt and depth", "(b) inward flow (away from the wall)", "(c) defect motion vs own jet", "(d) Local wall flow, along the circulation")):
    a.set_xlabel("time relative to departure (t.u.)"); a.set_ylabel(yl); a.set_title(ti, fontsize=10); a.axvline(0, color="0.6", lw=0.8, ls=":"); a.axhline(0, color="0.8", lw=0.6)
    a.legend(fontsize=8.5, frameon=False)
fig.tight_layout(w_pad=3); fig.savefig(os.path.join(H, "fig_compare_bulk.png"), dpi=200); print("ok")
