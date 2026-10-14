#!/usr/bin/env python3
"""Tilt law across lambda in active circles (outline plan A; runs/lambda_circle). Per run (second half): wall +1/2 count, bulk charge, wall-defect tilt tau (tail - outward normal,
wrapped to (-90, 90] deg for the nematic tail), median |tau|, tilt coherence |<sign tau>|, wall circulation and v_rms, single-vortex check via the sign of the wall circulation.
Grouped by (k, lambda); compared with arccos(1/lambda). Usage: tiltlaw.py RUNDIR [...] -> tiltlaw.txt, fig_tiltlaw.png"""
import os, sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
def one(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; lam = J["species"][0].get("tumble", 2.0); sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Tn /= np.abs(Tn); Nn = -1j * Tn
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn; Tn = -Tn
    t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5); h = len(fr) // 2; nw, qb, taus, coh = [], [], [], []
    for f in fr[h:]:
        if not len(f): continue
        d, j = tree.query(f[:, :2]); w = (np.sign(f[:, 2]) == sg) & (d < 4); nw.append(w.sum()); qb.append(0.5 * np.sign(f[d >= 4, 2]).sum())
        tau = np.angle(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))); taus += list(tau); coh.append(abs(np.mean(np.sign(tau))) if w.any() else np.nan)
    rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; ins = Path(W).contains_points(np.c_[X, Y]); dW, jW = tree.query(np.c_[X, Y]); vr, vt = [], []
    for tt in ft[len(ft) // 2:]:
        ii = idx[tt]; m = ins[ii]; vr.append(np.sqrt(np.mean(rows[ii, 4][m] ** 2 + rows[ii, 5][m] ** 2))); b = m & (dW[ii] > 1) & (dW[ii] < 5)
        vt.append(np.mean(rows[ii, 4][b] * Tn[jW[ii][b]].real + rows[ii, 5][b] * Tn[jW[ii][b]].imag))
    taus = np.degrees(np.abs(np.array(taus)))
    return (k, lam, np.mean(nw), np.mean(qb), np.median(taus), np.mean(taus), np.nanmean(coh), np.mean(np.abs(vt)), np.mean(vr), np.mean(np.sign(vt) == np.sign(np.mean(vt))), os.path.basename(run))
if __name__ == "__main__":
    runs = [os.path.abspath(r) for r in sys.argv[1:] if os.path.exists(r + "/DONE")]
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=1)
    L = ["k  lambda | runs | wall +1/2 (2k-2) | bulk charge | tilt |tau| median / mean (deg) | arccos(1/lambda) | tilt coherence | |wall v_t| | v_rms | P(circulation keeps its sign)"]
    G = {}
    for r in res: G.setdefault((r[0], r[1]), []).append(r)
    for key in sorted(G):
        S = np.array([r[2:10] for r in G[key]], float); m = S.mean(0); e = S.std(0) / np.sqrt(max(len(S) - 1, 1)); ac = np.degrees(np.arccos(1 / key[1])) if key[1] >= 1 else np.nan
        L.append(f"{key[0]:2d}  {key[1]:3.1f} | {len(S)} | {m[0]:5.2f}±{e[0]:.2f} ({2 * key[0] - 2}) | {m[1]:+.2f}±{e[1]:.2f} | {m[2]:5.1f} / {m[3]:5.1f} ± {e[3]:.1f} | {ac:5.1f} | {m[4]:.2f} | {m[5]:.3f} | {m[6]:.3f} | {m[7]:.2f}")
    print("\n".join(L)); open(os.path.join(H, os.environ.get("OUT", "tiltlaw.txt")), "w").write("\n".join(L) + "\n")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4)); lx = np.linspace(1.0, 6.5, 200)
    ax[0].plot(lx, np.degrees(np.arccos(1 / lx)), "k-", label="arccos(1/λ)")
    for kk, mk, col in ((4, "s", "#1f5fa8"), (6, "o", "#c2452d"), (8, "^", "#2a7a3b")):
        ks = sorted(l for (k_, l) in G if k_ == kk)
        if not ks: continue
        M = [np.array([r[2:10] for r in G[(kk, l)]], float) for l in ks]
        ax[0].errorbar(ks, [x[:, 3].mean() for x in M], [x[:, 3].std() / np.sqrt(max(len(x) - 1, 1)) for x in M], fmt=mk + "-", color=col, ms=6, capsize=3, label=f"MPCD k = {kk}")
        ax[1].plot(ks, [x[:, 0].mean() - (2 * kk - 2) for x in M], mk + "-", color=col, label=f"k = {kk}"); ax[2].plot(ks, [x[:, 4].mean() for x in M], mk + "-", color=col, label=f"k = {kk}")
    ax[0].set_xlabel("flow-alignment parameter λ"); ax[0].set_ylabel("wall-defect tilt |τ| (deg, mean)"); ax[0].legend(frameon=False, fontsize=8)
    ax[1].set_xlabel("λ"); ax[1].set_ylabel("wall +½ minus 2k−2"); ax[1].axhline(0, color="0.7", lw=0.8); ax[1].legend(frameon=False, fontsize=8)
    ax[2].set_xlabel("λ"); ax[2].set_ylabel("tilt coherence"); ax[2].legend(frameon=False, fontsize=8)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.suptitle("Active circles across λ (area 2500, natural activity): tilt law and the overcharged state", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(H, "fig_tiltlaw.png"), dpi=200)
