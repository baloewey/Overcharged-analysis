#!/usr/bin/env python3
"""Are the near-radial / wrong-sign wall defects of crowded layers failed escape run-ups? (2026-10-22) Shape-sweep events (shapes/escape_events/events.npz:
every wall +1/2 start, frames every 10 t.u.). Per (run, frame): majority tilt sign; classes by tilt relative to the majority: MAJ (> 20 deg), NEAR-RADIAL
(-20..20), WRONG (< -20). Per class and rho: share, mean depth, prompt escape probability (20 t.u.), lagged escape (in wall at +40, out at +60).
-> failed_runups.txt"""
import csv, numpy as np
S = "/home/baloewe/Physics/Overcharged/reports/master/theory/shapes"
E = np.load(S + "/escape_events/events.npz")["E"]; run = E[:, 0].astype(int); fr = E[:, 1].astype(int); tl = E[:, 2]; dep = E[:, 5]; esc = E[:, 7] > 0.5; late = E[:, 9] > 0.5
meta = list(csv.DictReader(open(S + "/escape_events/runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
rho = np.array([float(SR[m["run"]]["rho"]) for m in meta]); asp = np.array([float(SR[m["run"]]["aspect"]) for m in meta]); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
lev = np.argmin(np.abs(rho[:, None] - LV), 1)
key = run.astype(np.int64) * 100000 + fr; u, inv = np.unique(key, return_inverse=True); maj = np.sign(np.bincount(inv, np.sign(tl))); maj[maj == 0] = 1
rel = tl * maj[inv]; cls = np.where(rel > 20, 0, np.where(rel >= -20, 1, 2)); NM = ("majority", "near-radial", "wrong")
L = ["FAILED RUN-UPS? classes by tilt relative to the frame majority: majority > 20 deg, near-radial -20..20, wrong < -20. Per class: share | mean depth | p prompt | p lagged (x1000)"]
for nm, sel in (("circles", asp < 1.05), ("aspect>=1.4", asp >= 1.4)):
    L.append(f"\n{nm}")
    for i, r0 in enumerate(LV):
        m = np.isin(run, np.flatnonzero((lev == i) & sel)); row = []
        for c in range(3):
            q = m & (cls == c); row.append(f"{NM[c]}: {q.sum() / m.sum():.3f} | {dep[q].mean():.2f} | {1000 * esc[q].mean():5.1f} | {1000 * late[q].mean():5.1f}")
        L.append(f"  rho {r0:.2f}: " + "  ||  ".join(row))
open("failed_runups.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
