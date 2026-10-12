#!/usr/bin/env python3
"""How much of the inward flow at an escaping wall defect is its own? (user, 2026-10-11). Fine-flow runs (flow and defects every 2 t.u.), depleting escapes (energetics.events),
aligned on departure, lags -20..+2 frames. At each lag, for the escaper:
  depth; cos(tau); fluid inward velocity at the core (cells within 2) and in the ring 4-8 cells around it (inward = along -outward normal of the nearest wall point);
  the defect's own inward velocity (central difference of its depth along the chain, +-1 frame).
Self-propulsion calibration (same runs): bulk +1/2 (depth > 8) displacement over 2 t.u. along the head direction (-angle) minus the local fluid velocity at their core along the same
direction -> v0 (self-propulsion relative to the fluid). Prediction for the escaper's own inward contribution: v0 * cos(tau).
Usage: selfflow.py RUNDIR [...] -> selfflow.txt, fig_selfflow.png"""
import gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
LAGS = list(range(-20, 3))
def load_flow(run):
    fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"; op = gzip.open if fn.endswith(".gz") else open
    rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit())); ft = np.unique(rows[:, 0])
    idx = {tt: np.flatnonzero(rows[:, 0] == tt) for tt in ft}; return rows, ft, idx
def analyse(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=50); Q = bc["Q"]; sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]
    def frame_flow(tt):
        t0 = ft[np.argmin(np.abs(ft - tt))]; ii = idx[t0]; return X[ii], Y[ii], rows[ii, 4], rows[ii, 5]
    # v0 calibration from bulk +1/2
    v0s = []
    for n in range(0, len(fr) - 1, 5):
        f, g = fr[n], fr[n + 1]
        if not len(f) or not len(g): continue
        d, _ = tree.query(f[:, :2]); b = np.flatnonzero((f[:, 2] > 0) & (d > 8)); gp = g[g[:, 2] > 0, :2]
        if not len(b) or not len(gp): continue
        x, y, vx, vy = frame_flow(t[n])
        for i in b:
            dd = np.hypot(*(gp - f[i, :2]).T); kk = np.argmin(dd)
            if dd[kk] > 2: continue
            head = np.array([-np.cos(f[i, 3]), -np.sin(f[i, 3])]); disp = (gp[kk] - f[i, :2]) / (t[n + 1] - t[n])
            rr = np.hypot(x - f[i, 0], y - f[i, 1]); m = rr < 2; rg = (rr > 4) & (rr < 8)
            if not m.any() or not rg.any(): continue
            v0s.append(np.dot(np.array([vx[m].mean() - vx[rg].mean(), vy[m].mean() - vy[rg].mean()]), head))   # own jet: core minus ring flow, along the head
    out = []
    for e in ev:
        n = e["n"]; chain = dict(e["chain"]); c = chain[n]
        for m in range(n + 1, n + 4):
            h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
            if not len(hp): break
            dd = np.hypot(*(hp - c).T); kk = np.argmin(dd)
            if dd[kk] > 2.5: break
            c = hp[kk]; chain[m] = c
        depth = {m: tree.query(p)[0] for m, p in chain.items()}
        for lag in LAGS:
            m = n + lag
            if m not in chain: continue
            p = chain[m]; d, j = tree.query(p); nin = -Nn[j]; x, y, vx, vy = frame_flow(t[m]); r = np.hypot(x - p[0], y - p[1])
            core = r < 2; ring = (r > 4) & (r < 8)
            vin_core = np.mean(vx[core] * nin.real + vy[core] * nin.imag) if core.any() else np.nan
            vin_ring = np.mean(vx[ring] * nin.real + vy[ring] * nin.imag) if ring.any() else np.nan
            f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - p).T))); ct = np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j])))))
            vdef = (depth[m + 1] - depth[m - 1]) / (t[m + 1] - t[m - 1]) if (m + 1 in depth and m - 1 in depth) else np.nan
            out.append((lag, d, ct, vin_core, vin_ring, vdef))
    return np.array(out).reshape(-1, 6), np.array(v0s)
if __name__ == "__main__":
    O, V = [], []
    for run in sys.argv[1:]:
        try: o, v = analyse(os.path.abspath(run)); O.append(o); V.append(v); print(run, len(o), "rows,", len(v), "v0 samples", flush=True)
        except Exception as ex: print(run, ex)
    O = np.concatenate(O); V = np.concatenate(V); v0 = np.median(V); L = [f"own jet of bulk +1/2 (core minus ring flow along the head): median {v0:.3f}, mean {V.mean():.3f} cells/t.u. (n {len(V)})"]
    L.append("lag(t.u.)  depth  cos_tau  fluid_in_core  fluid_in_ring(4-8)  core-ring  defect_in_velocity  v0*cos_tau")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = []
    for lag in LAGS:
        a = O[O[:, 0] == lag]
        if len(a) < 5: continue
        mu = np.nanmean(a, 0); rows.append((2 * lag, mu[1], mu[2], mu[3], mu[4], mu[3] - mu[4], mu[5], v0 * mu[2]))
        L.append(f"{2 * lag:+5d}   {mu[1]:5.2f}  {mu[2]:+.3f}   {mu[3]:+.3f}        {mu[4]:+.3f}            {mu[3] - mu[4]:+.3f}     {mu[5]:+.3f}        {v0 * mu[2]:+.3f}")
    R = np.array(rows); plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(R[:, 0], R[:, 3], "o-", color="#c2452d", ms=3, label="fluid at the defect (within 2 cells)"); ax[0].plot(R[:, 0], R[:, 5], ":", color="#c2452d", label="at the defect minus ring"); ax[0].plot(R[:, 0], R[:, 4], "s-", color="#1f5fa8", ms=3, label="fluid in a ring 4-8 cells around it")
    ax[0].plot(R[:, 0], R[:, 6], "^-", color="0.3", ms=3, label="defect's own inward velocity"); ax[0].plot(R[:, 0], R[:, 7], "--", color="#3a8a4f", label=r"own-jet estimate $u_\mathrm{self}\cos\tau$ (bulk +1/2 calibration)")
    ax[0].axvline(0, color="0.6", lw=0.8, ls=":"); ax[0].set_xlabel("time relative to departure (t.u.)"); ax[0].set_ylabel("towards the cavity centre (cells/t.u.)"); ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title("(a) inward velocities at an escaping defect", fontsize=9)
    ax[1].plot(R[:, 0], R[:, 1], "o-", color="0.3", ms=3); ax[1].set_xlabel("time relative to departure (t.u.)"); ax[1].set_ylabel("depth (cells)"); ax1 = ax[1].twinx()
    ax1.plot(R[:, 0], R[:, 2], "s-", color="#d08a1e", ms=3); ax1.set_ylabel(r"$\cos\tau$", color="#d08a1e"); ax[1].set_title("(b) depth (black) and tilt (orange)", fontsize=9); ax[1].axvline(0, color="0.6", lw=0.8, ls=":")
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_selfflow{os.environ.get('TAG', '')}.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, f"selfflow{os.environ.get('TAG', '')}.txt"), "w").write("\n".join(L) + "\n")
