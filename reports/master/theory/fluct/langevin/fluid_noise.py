#!/usr/bin/env python3
"""Fluid velocity fluctuations felt by a wall defect (passive cavity, k = 6, R = 27.28; runs/fluidnoise/passive_R27_k06, coarseOut every 1 time unit).
Cell centre-of-mass velocities are smoothed with a Gaussian of width SIG = 1.5 cells (the defect-core scale, as the Stokes force smearing) and
sampled at 72 points on rings at depth d from the wall; decomposed into tangential and normal (outward) components.
Per depth: variance of each component (over points and frames) and its 1/e autocorrelation time. Writes fluid_noise.json."""
import gzip, json, os
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
H = os.path.dirname(os.path.abspath(__file__))
RUN = "/home/baloewe/Physics/Overcharged/runs/fluidnoise/passive_R27_k06"
bc = json.load(open(RUN + "/input.json"))["BC"][0]; Q = bc["Q"]; R = bc["R"]
fn = RUN + "/coarsegrain.dat.gz" if os.path.exists(RUN + "/coarsegrain.dat.gz") else RUN + "/coarsegrain.dat"
op = gzip.open if fn.endswith(".gz") else open
rows = np.loadtxt((l for l in op(fn, "rt") if l[:1].isdigit() and len(l.split()) >= 8), usecols=(0, 1, 2, 4, 5, 7))   # skip header lines
ts = np.unique(rows[:, 0]); nx = int(rows[:, 1].max()) + 1; ny = int(rows[:, 2].max()) + 1
nf = len(rows) // (nx * ny); rows = rows[:nf * nx * ny].reshape(nf, nx, ny, 6)
SIG = 1.5; depths = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]; ph = np.linspace(0, 2 * np.pi, 72, endpoint=False)
out = {"R": R, "frames": nf, "dt": float(ts[1] - ts[0]), "sigma_smooth": SIG, "depth": {}}
inside = (np.add.outer(np.arange(nx) + 0.5 - Q[0], np.zeros(ny)) ** 2 + np.add.outer(np.zeros(nx), np.arange(ny) + 0.5 - Q[1]) ** 2) < R ** 2
VT = {d: [] for d in depths}; VN = {d: [] for d in depths}
for f in range(nf):
    vx = rows[f, :, :, 3] * inside; vy = rows[f, :, :, 4] * inside; w = gaussian_filter(inside.astype(float), SIG)
    sx = gaussian_filter(vx, SIG) / np.maximum(w, 1e-9); sy = gaussian_filter(vy, SIG) / np.maximum(w, 1e-9)
    for d in depths:
        r = R - d; X = Q[0] + r * np.cos(ph) - 0.5; Y = Q[1] + r * np.sin(ph) - 0.5     # array index coordinates (cell i centred at i + 0.5)
        ux = map_coordinates(sx, [X, Y], order=1); uy = map_coordinates(sy, [X, Y], order=1)
        VT[d].append(-np.sin(ph) * ux + np.cos(ph) * uy); VN[d].append(np.cos(ph) * ux + np.sin(ph) * uy)
def tau(A, dt):
    A = A - A.mean(0); c0 = (A * A).mean()
    for L in range(1, len(A) // 2):
        c = (A[:-L] * A[L:]).mean() / c0
        if c < np.exp(-1):
            cp = (A[:-L + 1] * A[L - 1:]).mean() / c0 if L > 1 else 1.0
            return dt * (L - 1 + (cp - np.exp(-1)) / (cp - c))
    return np.nan
for d in depths:
    T = np.array(VT[d]); Nn = np.array(VN[d])
    out["depth"][str(d)] = dict(var_t=float(T.var()), var_n=float(Nn.var()), tau_t=float(tau(T, out["dt"])), tau_n=float(tau(Nn, out["dt"])))
    print(f"depth {d}: rms tangential {np.sqrt(T.var()):.4f}, normal {np.sqrt(Nn.var()):.4f}; tau_t {out['depth'][str(d)]['tau_t']:.2f}, tau_n {out['depth'][str(d)]['tau_n']:.2f}")
json.dump(out, open(os.path.join(H, "fluid_noise.json"), "w"), indent=1)
