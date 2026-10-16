#!/usr/bin/env python3
"""Why do crowded wall layers tilt less? (user, 2026-10-15; RESEARCH.md 'TILT FALLS IN CROWDED LAYERS'). Main circle sweep, act 0.07, lambda 2, +k (2-10),
0.8 <= R/l_a <= 2.2, all 10 seeds (runs/cavity/server_main/main_act0.07). Defect frames every 2 t.u., first 20% dropped.
Per frame, wall +1/2 (within 4 cells of the wall): tilt tau_j = psi_j - phi_j (psi = defect angle, phi = polar angle = outward normal), wrapped to (-pi, pi];
frame coherence C = |<exp(i tau)>|, frame mean tilt <tau> = arg<exp(i tau)>, its sign s = sign(<tau>) (circulation sense, without the flow).
Per run: tilt of the coherent layer = median over frames of |<tau>| among frames with C > 0.9 (a measurement where the layer clearly shares one tilt);
fraction of frames with C > 0.9; per-defect tilt s*tau_j against the arc distance to the nearest wall neighbour (in units of the mean spacing 2 pi R/(2k-2)),
pooled per rho bin -> does an individual defect with a close neighbour tilt less, at fixed rho?
-> tilt_crowding.txt, fig_tilt_crowding.png (+pdf). Usage: tilt_crowding.py [NPROC]"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; n0 = 2 * k - 2; a = 2 * np.pi / n0
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2)
    Cs, Ms, NN = [], [], []
    for f in fr:
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        if w.sum() < 3: continue
        phi = np.arctan2(f[w, 1], f[w, 0]); tau = np.angle(np.exp(1j * (f[w, 3] - phi))); z = np.mean(np.exp(1j * tau))
        Cs.append(abs(z)); Ms.append(np.angle(z)); s = np.sign(np.angle(z))
        o = np.argsort(phi); ph = phi[o]; dphi = np.diff(np.r_[ph, ph[0] + 2 * np.pi]); nn = np.minimum(dphi, np.roll(dphi, 1)) / a   # nearest-neighbour arc / mean spacing
        NN += list(zip(nn, s * tau[o]))
    Cs, Ms = np.array(Cs), np.array(Ms); coh = Cs > 0.9
    return dict(run=run, k=k, R=R, rho=k / (np.pi * R), C=np.median(Cs), fcoh=coh.mean(), tilt_coh=np.degrees(np.median(np.abs(Ms[coh]))) if coh.sum() > 20 else np.nan,
                tilt_all=np.degrees(np.median(np.abs(Ms))), nn=np.array(NN, dtype=np.float32))


if __name__ == "__main__":
    runs = []
    for d in sorted(glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/DONE")):
        d = os.path.dirname(d); bc = json.load(open(d + "/input.json"))["BC"][0]
        if 2 <= bc["customAnch"] <= 10 and 0.8 <= bc["R"] / LA <= 2.2: runs.append(d)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs, chunksize=4)
    np.savez(os.path.join(H, "tilt_crowding_runs.npz"), **{c: np.array([x[c] for x in res]) for c in ("k", "R", "rho", "C", "fcoh", "tilt_coh", "tilt_all")},
             run=np.array([os.path.relpath(x["run"], ROOT) for x in res]))
    rb = [(0.025, 0.055), (0.055, 0.075), (0.075, 0.095), (0.095, 0.11), (0.11, 0.13), (0.13, 0.16)]
    L = [f"TILT vs WALL DENSITY (circles, lambda 2, k 2-10, 0.8 <= R/l_a <= 2.2, {len(res)} runs, 10 seeds). Tilt = |frame mean of tau| (deg); C = frame coherence |<exp(i tau)>|.",
         "rho bin | runs | median C | frac. frames C > 0.9 | tilt, all frames | tilt, coherent frames (C > 0.9) | per-defect signed tilt by nearest-neighbour arc (/mean spacing): <0.6 / 0.6-0.9 / 0.9-1.1 / >1.1 (deg; share of defects)"]
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2)); cols = plt.cm.viridis(np.linspace(0, 0.9, len(rb)))
    for (lo, hi), col in zip(rb, cols):
        S = [x for x in res if lo <= x["rho"] < hi]
        if not S: continue
        nn = np.concatenate([x["nn"] for x in S]); bins = [(0, 0.6), (0.6, 0.9), (0.9, 1.1), (1.1, 9)]
        per = [(np.degrees(np.median(nn[(nn[:, 0] >= a) & (nn[:, 0] < b), 1])), np.mean((nn[:, 0] >= a) & (nn[:, 0] < b))) for a, b in bins]
        tc = np.nanmedian([x["tilt_coh"] for x in S])
        L.append(f"[{lo:.3f},{hi:.3f}) | {len(S):3d} | {np.median([x['C'] for x in S]):.2f} | {np.median([x['fcoh'] for x in S]):.2f} | {np.median([x['tilt_all'] for x in S]):5.1f} | "
                 f"{tc:5.1f} (runs with >20 coherent frames: {np.sum(np.isfinite([x['tilt_coh'] for x in S]))}) | " + " / ".join(f"{p:5.1f} ({sh:.2f})" for p, sh in per))
        edges = np.linspace(0.2, 1.6, 15); mid = 0.5 * (edges[1:] + edges[:-1])
        med = [np.degrees(np.median(nn[(nn[:, 0] >= a) & (nn[:, 0] < b), 1])) if ((nn[:, 0] >= a) & (nn[:, 0] < b)).sum() > 200 else np.nan for a, b in zip(edges[:-1], edges[1:])]
        ax[2].plot(mid, med, "o-", ms=3, color=col, label=f"$\\rho$ {lo:.3f}-{hi:.3f}")
    rho = np.array([x["rho"] for x in res])
    ax[0].scatter(rho, [x["tilt_all"] for x in res], s=6, c="0.6", label="all frames")
    ax[0].scatter(rho, [x["tilt_coh"] for x in res], s=6, c="#c2452d", label="coherent frames ($C>0.9$)")
    ax[0].axhline(60, color="k", ls="--", lw=0.8); ax[0].set_xlabel("wall density $\\rho=k/\\pi R$"); ax[0].set_ylabel("tilt of the wall layer (deg)"); ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title("(a) Tilt vs wall density", fontsize=10)
    ax[1].scatter(rho, [x["C"] for x in res], s=6, c="0.3"); ax[1].set_xlabel("wall density $\\rho=k/\\pi R$"); ax[1].set_ylabel("median frame coherence $C$"); ax[1].set_title("(b) Coherence", fontsize=10)
    ax[2].axvline(1, color="0.8", lw=0.8); ax[2].set_xlabel("nearest-neighbour arc / mean spacing"); ax[2].set_ylabel("signed tilt of a wall defect (deg, median)")
    ax[2].set_title("(c) Does a close neighbour reduce the tilt?", fontsize=10); ax[2].legend(fontsize=7, frameon=False)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_tilt_crowding.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "tilt_crowding.txt"), "w").write("\n".join(L) + "\n")
