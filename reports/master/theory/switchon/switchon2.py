#!/usr/bin/env python3
"""Passive -> active switch across crowding and shape (runs/switchon2; 2026-10-22). Circles R 31.83, k 4 / 10 / 16 (rho 0.04 / 0.10 / 0.16); ellipses aspect
1.5 and 2 (area 2500, k 6); initial director along x or random. Passive 3000 t.u. (after warm-up), then active (zeta 0.07) from the checkpoint.
Per run, 20 t.u. bins over both stages: wall +1/2 count (d < 4), tilt tau = circular mean of (tail - outward normal) over wall +1/2, coherence |<e^{i tau}>|,
circles: registry arg Z_N (N = 2k-2, Z_N = <exp(i N phi)>) and lock offset D = arg Z_N - tau over frames with exactly N wall defects; circulation <v_t>/<|v|>
about the centroid from the flow (active stage). Final state = last 1000 t.u. of the active stage; switch-on time = first bin after the switch where the
coherence exceeds 0.9 (circles) and |tau| reaches 80% of its final value. -> switchon2.txt, switchon2.npz"""
import glob, gzip, json, os, sys, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle; from wall_fourier import wall_polyline
RD = "/home/baloewe/Physics/Overcharged/runs/switchon2"; BIN = 20.0
def flow_circ(run, c):
    T, C = [], []; cur, num, den = None, 0.0, 0.0
    with gzip.open(run + "/flowfield.dat.gz", "rt") as f:
        for line in f:
            p = line.split()
            if len(p) != 7 or not p[0][0].isdigit(): continue
            t = float(p[0]); x, y = float(p[1]) + 0.5 - c[0], float(p[2]) + 0.5 - c[1]; vx, vy = float(p[4]), float(p[5])
            if cur is not None and t != cur: T.append(cur); C.append(num / max(den, 1e-12)); num = den = 0.0
            cur = t; r = np.hypot(x, y)
            if r > 0.5 and (vx or vy): num += (x * vy - y * vx) / r; den += np.hypot(vx, vy)
    if cur is not None: T.append(cur); C.append(num / max(den, 1e-12))
    return np.array(T), np.array(C)
def stage(run, bc, W, tree, Nn, N, circ):
    c = np.array(bc["Q"][:2], float); t, fr = read_with_angle(run, c, 0.0); rows = []
    for tt, f in zip(t, fr):
        if not len(f): continue
        p = f[f[:, 2] > 0]; d, j = tree.query(p[:, :2]); w = d < 4
        if not w.any(): rows.append((tt, 0, np.nan, np.nan, np.nan, np.nan)); continue
        tl = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]])))); z = np.mean(np.exp(1j * tl)); ph = np.arctan2(p[w, 1], p[w, 0])
        zn = np.mean(np.exp(1j * N * ph)) if circ else np.nan
        rows.append((tt, w.sum(), np.angle(z), abs(z), np.angle(zn) if circ else np.nan, np.angle(zn * np.conj(z)) if circ and w.sum() == N else np.nan))
    return np.array(rows)
out = []; L = ["SWITCH-ON SCAN (passive -> active, zeta 0.07): per run, passive end (last 500 t.u.) and active final (last 1000 t.u.)",
               "run | k | N=2k-2 | rho | passive: wall, C, |tau| | active final: wall, C, |tau| (deg), lock offset D (circ. mean, deg; spread) | circulation | switch-on time (t.u.)"]
for name in sorted(os.listdir(RD)):
    if not (os.path.isdir(f"{RD}/{name}/active") and os.path.exists(f"{RD}/{name}/passive/input.json")): continue
    bc = json.load(open(f"{RD}/{name}/passive/input.json"))["BC"][0]; k = bc["customAnch"]; N = 2 * k - 2; circ = name.startswith("circle")
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; T_ = np.roll(Wc, -1) - np.roll(Wc, 1); T_ /= np.abs(T_); Nn = -1j * T_
    if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
    A = stage(f"{RD}/{name}/passive", bc, W, tree, Nn, N, circ); B = stage(f"{RD}/{name}/active", bc, W, tree, Nn, N, circ)
    ts = B[0, 0]; pa = A[A[:, 0] > A[-1, 0] - 500]; fa = B[B[:, 0] > B[-1, 0] - 1000]
    tf, cf = flow_circ(f"{RD}/{name}/active", np.array(bc["Q"][:2], float)); circv = np.mean(np.abs(cf[tf > tf[-1] - 1000])) if len(cf) else np.nan
    tauF = np.degrees(np.nanmedian(np.abs(fa[:, 2]))); Dm = np.angle(np.nanmean(np.exp(1j * fa[:, 5]))) if circ else np.nan
    Ds = np.degrees(np.sqrt(-2 * np.log(max(abs(np.nanmean(np.exp(1j * fa[:, 5]))), 1e-9)))) if circ and np.isfinite(fa[:, 5]).any() else np.nan
    sw = np.nan
    for i in range(0, len(B) - 10, 10):
        seg = B[i:i + 10]
        if np.nanmean(seg[:, 3]) > (0.9 if circ else 0.6) and np.nanmedian(np.abs(np.degrees(seg[:, 2]))) > 0.8 * tauF: sw = seg[0, 0] - ts; break
    L.append(f"{name:16s} | {k:2d} | {N:2d} | {2 * k / P:.3f} | {np.nanmean(pa[:, 1]):5.1f}, {np.nanmean(pa[:, 3]):.2f}, {np.degrees(np.nanmedian(np.abs(pa[:, 2]))):5.1f} | "
             f"{np.nanmean(fa[:, 1]):5.1f}, {np.nanmean(fa[:, 3]):.2f}, {tauF:5.1f}, D {np.degrees(Dm):+6.1f} (spread {Ds:5.1f}) | {circv:.2f} | {sw:6.0f}")
    out.append((name, A, B))
L.append("\nRESULTS (filled in after inspection below).")
open("switchon2.txt", "w").write("\n".join(L) + "\n"); np.savez("switchon2.npz", **{f"{n}_passive": a for n, a, b in out}, **{f"{n}_active": b for n, a, b in out}); print("\n".join(L))
