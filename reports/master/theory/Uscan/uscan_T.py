#!/usr/bin/env python3
"""Near-wall temperature in the U-scan (2026-10-27). As crowd_protect/drive_barrier/analyse_db.py (UMODE relaxed): per run, fit ln S(x) = a - U_K(x)/T over
x >= d0 + 0.5 with S n >= 5, where U_K = (K_U / 44) x the relaxed-neighbour one-constant drag energy (drag_relaxed.npz, K 44, N = rho x 200).
K_U from the fluctuation splay constant (frank_U.txt, dir_ml, qmax 0.5, mean of 2 seeds): U 5 -> 8.4, U 10 -> 43.6, U 20 -> 89.3 kBT.
If the barrier sets the escape and the drive does not depend on K, T should be the same at every U. -> uscan_T.txt"""
import csv, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); Z = np.load(os.path.join(H, "..", "crowd_protect", "drag_relaxed.npz"))
NS = np.array([8, 12, 18, 24, 28]); DEL = Z["DEL"]; EN = np.array([Z[f"relaxed_N{N}"] for N in NS]); KU = {5: 8.4, 10: 43.6, 20: 89.3}
def U(x, d0, rho):
    N = np.clip(rho * 200, NS[0], NS[-1]); E = np.array([np.interp(N, NS, EN[:, j]) for j in range(len(DEL))]); return np.interp(np.clip(x - d0, 0, DEL[-1]), DEL, E)
R = [r for r in csv.DictReader(open(os.path.join(H, "uscan_escape.csv"))) if not r.get("err")]
XG = np.array(sorted(float(k[1:]) for k in R[0] if k.startswith("S") and k[1:2].isdigit()))
out = {}
for r in R:
    Uv, k, n, d0 = int(r["U"]), int(r["k"]), int(r["n"]), float(r["d0"]); rho = {5: 0.04, 14: 0.12}[k]
    S = np.array([float(r[f"S{x:g}"]) for x in XG]); ok = (XG >= d0 + 0.5) & (S * n >= 5) & (S > 0)
    u = np.array([U(x, d0, rho) for x in XG[ok]]) * KU[Uv] / 44.0; y = np.log(S[ok]); b = np.polyfit(u, y, 1)
    out.setdefault((k, Uv), []).append((-1 / b[0] if b[0] < 0 else np.nan, 1 - np.var(y - np.polyval(b, u)) / max(np.var(y), 1e-12)))
L = ["U-SCAN NEAR-WALL TEMPERATURE (aspect-2 ellipse, R_eq 1.5 l_a; U_K = K/44 x relaxed one-constant drag energy). Mean +- s.e. over seeds.",
     "k (rho) | mfpot | K (kBT) | seeds | T (kBT) | T/K | R^2"]
for (k, Uv), v in sorted(out.items()):
    a = np.array(v); t = a[:, 0]
    L.append(f"{k:2d} ({ {5: 0.04, 14: 0.12}[k]:.2f}) | {Uv:2d} | {KU[Uv]:5.1f} | {len(a):2d} | {np.nanmean(t):5.1f} +- {np.nanstd(t, ddof=1) / np.sqrt(len(t)):.1f} | {np.nanmean(t) / KU[Uv]:.3f} | {np.nanmean(a[:, 1]):.2f}")
open(os.path.join(H, "uscan_T.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
