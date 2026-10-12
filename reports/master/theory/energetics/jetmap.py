#!/usr/bin/env python3
"""Own-jet strength of +1/2 defects vs depth (distance to the wall) and tilt (user, 2026-10-11). All +1/2 in sampled frames (every 3rd flow frame, second half), fine-flow runs.
jet vector = mean flow within 2 cells minus mean flow in the 4-8 cell ring (ring restricted to fluid cells inside the cavity); projected (i) on the head (-angle) -> jet strength,
(ii) on the inward wall normal -> inward jet. Tilt: cos tau = cos(tail - outward normal) (1 = head into the bulk). Binned by depth (0-1, 1-2, ..., 8-10, >10) x cos tau (5 bins).
Usage: jetmap.py TAG RUNDIR [...] -> jetmap_<TAG>.txt, jetmap_<TAG>.npy, fig_jetmap_<TAG>.png"""
import os, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
from edge_criteria import read_with_angle  # noqa
DB = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 99.0]); CB = np.linspace(-1, 1, 6)
def one(run):
    J = __import__("json").load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t)
    rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; ins = Path(W).contains_points(np.c_[X, Y]); out = []
    for tt in ft[len(ft) // 2::3]:
        n = int(np.argmin(np.abs(t - tt)))
        if abs(t[n] - tt) > 0.6 or not len(fr[n]): continue
        ii = idx[tt]; x, y, vx, vy, ok = X[ii], Y[ii], rows[ii, 4], rows[ii, 5], ins[ii]
        f = fr[n]
        for d_, j, fd in zip(*tree.query(f[:, :2]), f):
            if fd[2] <= 0: continue
            r = np.hypot(x - fd[0], y - fd[1]); c = (r < 2) & ok; g = (r > 4) & (r < 8) & ok
            if not c.any() or g.sum() < 5: continue
            jv = np.array([vx[c].mean() - vx[g].mean(), vy[c].mean() - vy[g].mean()]); hd = -np.exp(1j * fd[3]); nin = -Nn[j]
            if not Path(W).contains_point(fd[:2]): d_ = -d_
            out.append((d_, np.cos(np.angle(np.exp(1j * (fd[3] - np.angle(Nn[j]))))), jv[0] * hd.real + jv[1] * hd.imag, jv[0] * nin.real + jv[1] * nin.imag))
    return out
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    if os.environ.get("REPLOT"): A = np.load(os.path.join(H, f"jetmap_{tag}.npy"))
    else:
        with Pool(min(12, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
        A = np.array(sum(res, [])); np.save(os.path.join(H, f"jetmap_{tag}.npy"), A)
    d, ct, jh, ji = A.T
    L = [f"{len(A)} +1/2 samples. Jet along the head (cells/t.u.), mean ± s.e. [n]; columns: depth bins; rows: cos tau bins (1 = head into the bulk)",
         "depth:         " + "  ".join(f"{DB[i]:g}-{DB[i+1]:g}".rjust(16) for i in range(len(DB) - 1))]
    M = np.full((len(CB) - 1, len(DB) - 1), np.nan); Mi = M.copy()
    for a in range(len(CB) - 1):
        cells = []
        for b in range(len(DB) - 1):
            m = (ct >= CB[a]) & (ct < CB[a + 1] + (1e-9 if a == len(CB) - 2 else 0)) & (d >= DB[b]) & (d < DB[b + 1])
            if m.sum() >= 20: M[a, b] = jh[m].mean(); Mi[a, b] = ji[m].mean()
            cells.append(f"{jh[m].mean():+.3f}±{jh[m].std() / np.sqrt(max(m.sum(), 1)):.3f}[{m.sum()}]".rjust(16) if m.sum() >= 20 else "-".rjust(16))
        L.append(f"cos {CB[a]:+.1f}..{CB[a+1]:+.1f}: " + "  ".join(cells))
    L.append("\nall tilts, by depth: jet along head | inward jet | n")
    for b in range(len(DB) - 1):
        m = (d >= DB[b]) & (d < DB[b + 1])
        if m.sum(): L.append(f"  depth {DB[b]:g}-{DB[b+1]:g}: {jh[m].mean():+.3f}±{jh[m].std() / np.sqrt(m.sum()):.3f} | {ji[m].mean():+.3f} | {m.sum()}")
    bulk = jh[d > 10].mean(); dd = np.array([0.5 * (DB[b] + min(DB[b + 1], 14)) for b in range(len(DB) - 1)])
    prof = np.array([jh[(d >= DB[b]) & (d < DB[b + 1])].mean() for b in range(len(DB) - 1)]) / bulk
    L.append("\njet strength j(d)/j_bulk: " + "  ".join(f"{x:g}:{y:.2f}" for x, y in zip(dd, prof)))
    gg = []
    for b in range(len(DB) - 1):
        m = (d >= DB[b]) & (d < DB[b + 1]); x = jh[m] * ct[m]; gg.append(np.sum(ji[m] * x) / np.sum(x * x) if m.sum() else np.nan)
    gg = np.array(gg); L.append("wall-normal transmission g(d) = <j_in * j_head cos tau> / <(j_head cos tau)^2>: " + "  ".join(f"{x:g}:{y:.2f}" for x, y in zip(dd, gg)))
    print("\n".join(L)); open(os.path.join(H, f"jetmap_{tag}.txt"), "w").write("\n".join(L) + "\n")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    cmap = plt.get_cmap("viridis")
    for a in range(len(CB) - 1):
        ax[0].plot(dd, M[a], "-o", ms=3, color=cmap(a / (len(CB) - 2)), label=f"cos τ {CB[a]:+.1f}…{CB[a+1]:+.1f}")
        ax[1].plot(dd, Mi[a], "-o", ms=3, color=cmap(a / (len(CB) - 2)))
    ax[0].set_ylabel("jet along the head (cells/t.u.)"); ax[1].set_ylabel("jet towards the cavity centre (cells/t.u.)")
    ax[2].plot(dd, gg, "-o", color="k", label="normal component transmitted, g"); ax[2].plot(dd, prof, "--s", color="0.5", ms=3, label="jet strength / bulk value")
    ax[2].set_ylabel("ratio"); ax[2].set_ylim(0, 1.4); ax[2].axhline(1, color="0.8", lw=0.6); ax[2].legend(fontsize=7, frameon=False)
    ax[0].set_title("(a) jet strength (along the head)", fontsize=9); ax[1].set_title("(b) jet towards the cavity centre", fontsize=9); ax[2].set_title("(c) wall blocks the normal component", fontsize=9)
    for a_ in ax: a_.set_xlabel("depth (cells from the wall)"); a_.axhline(0, color="0.8", lw=0.6); [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax[0].legend(fontsize=7, frameon=False); fig.suptitle(f"Own jet of +1/2 vs depth and tilt ({tag}; {len(A)} samples)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_jetmap_{tag}.png"), dpi=200)
