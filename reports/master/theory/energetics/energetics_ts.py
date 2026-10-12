#!/usr/bin/env python3
"""Time series around depleting escapes, aligned on departure (frame n): lags -30..+2 frames (2 t.u. each). For the escaper (followed backwards along its chain; forwards
for lags > 0) and for a random resident control in the same frame and region: elastic force F0, cost E(2)-E(0) and E(6)-E(0) (measured configuration, energetics.landscape),
measured flow at the probe (inward, |tangential|; nearest flow frame), and cos(tau). Usage: energetics_ts.py RUNDIR [...] -> ets_<tag>.npy (rows) and fig_ets_<tag>.png"""
import os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
LAGS = (-30, -25, -20, -15, -10, -7, -5, -3, -2, -1, 0, 1, 2)
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=10); eps = bc["aInv"][0] if "aInv" in bc else 1.0; E = En.bd.Ellipse(float(eps)); Q = bc["Q"]; sg = 1 if k > 0 else -1
    from scipy.spatial import cKDTree
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max(); rng = np.random.default_rng(11); out = []
    def probe_cos(frame, p):
        zs = frame[:, 0] + 1j * frame[:, 1]; i = int(np.argmin(np.abs(zs - (p[0] + 1j * p[1])))); nrm = E.Nn[np.argmin(np.abs(E.W - zs[i]))]
        return np.cos(np.angle(np.exp(1j * (frame[i, 3] - np.angle(nrm)))))
    for e in ev:
        n = e["n"]; chain = dict(e["chain"]); c = chain[n]
        for m in range(n + 1, n + 3):                                   # forwards
            h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
            if not len(hp): break
            d2 = np.hypot(*(hp - c).T); kk = np.argmin(d2)
            if d2[kk] > 2.5: break
            c = hp[kk]; chain[m] = c
        for lag in LAGS:
            m = n + lag
            if m not in chain: continue
            f = fr[m]; L = En.landscape(E, f, chain[m], k)
            if L is None: continue
            vx, vy = En.flow_at(run, Q, t[m], chain[m]); zp = chain[m][0] + 1j * chain[m][1]; nrm = E.Nn[np.argmin(np.abs(E.W - zp))]
            out.append(("esc", lag, L["F0"], L["prof"][3], L["E6"], -(vx * nrm.real + vy * nrm.imag), abs(-vx * nrm.imag + vy * nrm.real), probe_cos(f, chain[m])))
            d, j = tree.query(f[:, :2]); res = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)) if np.hypot(*(f[jj, :2] - chain[m])) > 3 and abs(abs(f[jj, 0]) - abs(chain[m][0])) / A < 0.15]
            if res and lag in (-30, -20, -10, -5, 0):
                pc = f[rng.choice(res), :2]; Lc = En.landscape(E, f, pc, k)
                if Lc is not None:
                    vx, vy = En.flow_at(run, Q, t[m], pc); zc = pc[0] + 1j * pc[1]; nc = E.Nn[np.argmin(np.abs(E.W - zc))]
                    out.append(("ctl", lag, Lc["F0"], Lc["prof"][3], Lc["E6"], -(vx * nc.real + vy * nc.imag), abs(-vx * nc.imag + vy * nc.real), probe_cos(f, pc)))
    return out
def plot(rows, tag, title):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    R = np.array([r[1:] for r in rows if r[0] == "esc"], float); C = np.array([r[1:] for r in rows if r[0] == "ctl"], float)
    panels = [(1, "elastic force on the defect\n(kBT/cell; < 0 = back to the wall)"), (2, "elastic cost to move 2 cells in (kBT)"), (3, "elastic cost to move 6 cells in (kBT)"),
              (4, "flow at the defect, inward (cells/t.u.)"), (5, "flow at the defect, along the wall (cells/t.u.)"), (6, r"tilt: $\cos\tau$ (1 = head into the bulk)")]
    fig, axs = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True); axs = axs.ravel()
    for ax, (q, lab) in zip(axs, panels):
        lags = sorted(set(R[:, 0])); mu = [np.nanmean(R[R[:, 0] == l, q]) for l in lags]; se = [np.nanstd(R[R[:, 0] == l, q]) / np.sqrt(max(np.isfinite(R[R[:, 0] == l, q]).sum(), 1)) for l in lags]
        x = 2 * np.array(lags); ax.fill_between(x, np.array(mu) - se, np.array(mu) + se, color="#c2452d", alpha=0.2, lw=0); ax.plot(x, mu, "o-", color="#c2452d", ms=3, lw=1.5, label="escaping defect")
        if len(C):
            cl = sorted(set(C[:, 0])); cm = [np.nanmean(C[C[:, 0] == l, q]) for l in cl]; cs = [np.nanstd(C[C[:, 0] == l, q]) / np.sqrt(max(np.isfinite(C[C[:, 0] == l, q]).sum(), 1)) for l in cl]
            ax.errorbar(2 * np.array(cl), cm, yerr=cs, fmt="s--", color="0.45", ms=3, lw=1, capsize=2, label="resident wall defects (same frames)")
        ax.axvline(0, color="0.6", lw=0.8, ls=":"); ax.set_title(lab, fontsize=8.5)
    for ax in axs[3:]: ax.set_xlabel("time relative to departure (t.u.)")
    axs[0].legend(fontsize=7, frameon=False); fig.suptitle(title, fontsize=9.5); fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_ets_{tag}.png"), dpi=200)
if __name__ == "__main__":
    runs = [os.path.abspath(r) for r in sys.argv[1:-1]]; tag = sys.argv[-1]
    with Pool(min(8, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    rows = sum(res, []); np.save(os.path.join(H, f"ets_{tag}.npy"), np.array(rows, dtype=object), allow_pickle=True)
    ne = len(set((r[1]) for r in rows)); plot(rows, tag, f"Depleting escapes, aligned on departure ({tag}; {sum(1 for r in rows if r[0] == 'esc' and r[1] == 0)} events)")
    print("done", tag, len(rows))
