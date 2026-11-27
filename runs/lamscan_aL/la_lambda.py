#!/usr/bin/env python3
"""Active length vs lambda (runs/lamscan_aL, L200 periodic, zeta 0.07, chi 0.5): l_a = sqrt(A / <N_defects>), N = all defects per frame, t >= half the run
(stationarity: first vs second quarter of that half). Run on the server in runs/lamscan_aL. -> la_lambda.txt"""
import gzip, glob, json, numpy as np
rows = []
for r in sorted(glob.glob("aL_lam*_s*")):
    J = json.load(open(r + "/input.json")); lam = J["species"][0]["tumble"]; A = J["domain"][0] * J["domain"][1]
    t, n = [], []
    for l in gzip.open(r + "/defects.dat.gz", "rt"):
        p = l.split()
        if len(p) == 2:
            try: t.append(float(p[0])); n.append(int(p[1]))
            except ValueError: pass
    t, n = np.array(t), np.array(n, float); m = t >= t[0] + 0.5 * (t[-1] - t[0]); nn = n[m]; h = len(nn) // 2
    rows.append((lam, r, nn.mean(), np.sqrt(A / nn.mean()), nn[:h].mean(), nn[h:].mean(), t[-1]))
L = ["lambda | run | <N> | l_a = sqrt(A/<N>) | <N> first / second quarter | t_end"]
for x in rows: L.append(f"{x[0]:3.1f} | {x[1]} | {x[2]:7.1f} | {x[3]:5.1f} | {x[4]:.1f} / {x[5]:.1f} | {x[6]:.0f}")
L.append("\nper lambda: l_a mean (seeds)")
for lam in sorted({x[0] for x in rows}):
    v = [x[3] for x in rows if x[0] == lam]; L.append(f"{lam:3.1f}: {np.mean(v):5.1f}  (" + ", ".join(f"{u:.1f}" for u in v) + ")")
open("la_lambda.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
