#!/usr/bin/env python3
"""Q criterion around defects (user, 2026-10-11; Head, Dore, ..., Shendruk, Nat. Phys. 20, 492 (2024): +1/2 defects sit on Q = 0 'viscometric' lines, P decays exponentially away
from them; -1/2 unconstrained). Q = 1/2((tr L)^2 - tr L^2) = det(grad u) in 2D (= 1/2(|Omega|^2 - |E|^2) if incompressible); normalised Qn = (|Omega|^2 - |E|^2)/(|Omega|^2 + |E|^2)
(-1 pure strain, 0 simple shear, +1 pure rotation), |Omega|^2 = w^2/2, |E|^2 = exx^2 + eyy^2 + 2 exy^2. Flow smoothed with a Gaussian of SIGMA cells (fluid cells inside the
cavity). Q = 0 set = cells whose Q sign differs from a 4-neighbour inside the cavity; d_Q = Euclidean distance to that set. Rows: depleting escapers (lags -30..+3 frames),
residents (lag % 5 == 0), all +1/2 and -1/2 (every 9th frame, second half), random fluid points (null, same frames). Usage: qcrit.py TAG RUNDIR [...] [env SIGMA=1]
-> qcrit_<TAG>.txt, qcrit_<TAG>.npy, fig_qcrit_<TAG>.png"""
import os, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, distance_transform_edt, map_coordinates
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
SIG = float(os.environ.get("SIGMA", 1.0)); LAGS = list(range(-30, 4))
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=50); Q = bc["Q"]; sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int); nx, ny = ix.max() + 1, iy.max() + 1
    XX, YY = np.meshgrid(np.arange(nx) + 0.5 - Q[0], np.arange(ny) + 0.5 - Q[1], indexing="ij"); inside = Path(W).contains_points(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape)
    dW = tree.query(np.c_[XX.ravel(), YY.ravel()])[0].reshape(XX.shape); cache = {}
    def field(tt):
        t0 = ft[np.argmin(np.abs(ft - tt))]
        if t0 in cache: return cache[t0]
        ii = idx[t0]; U = np.zeros((nx, ny)); V = np.zeros((nx, ny)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
        w8 = gaussian_filter(inside.astype(float), SIG) + 1e-9; Us = gaussian_filter(np.where(inside, U, 0), SIG) / w8; Vs = gaussian_filter(np.where(inside, V, 0), SIG) / w8
        ux, uy = np.gradient(Us); vx, vy = np.gradient(Vs); Qd = ux * vy - uy * vx; om = vx - uy; exy = 0.5 * (uy + vx)
        O2 = om ** 2 / 2; E2 = ux ** 2 + vy ** 2 + 2 * exy ** 2; Qn = (O2 - E2) / (O2 + E2 + 1e-12)
        sgn = np.sign(Qd); zero = np.zeros_like(inside)
        for sh in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nb = np.roll(sgn, sh, axis=(0, 1)); nbi = np.roll(inside, sh, axis=(0, 1)); zero |= inside & nbi & (nb != sgn)
        dQ = distance_transform_edt(~zero); Qd[~inside] = np.nan; Qn[~inside] = np.nan
        cache.clear(); cache[t0] = (Qd, Qn, dQ); return cache[t0]
    def at(F, p):   # bilinear interpolation at a point (cavity-centred coords -> lattice index)
        return float(map_coordinates(F, [[p[0] + Q[0] - 0.5], [p[1] + Q[1] - 0.5]], order=1, mode="nearest")[0])
    def samp(kind, q, lag, m, p, extra=np.nan):
        Qd, Qn, dQ = field(t[m]); d = tree.query(p)[0]
        if not Path(W).contains_point(p): return
        out.append((kind, q, lag, d, extra, at(np.nan_to_num(Qn), p), at(np.nan_to_num(Qd), p), at(dQ, p)))
    out = []; rng = np.random.default_rng(21); A = np.abs(W[:, 0]).max()
    def ctv(p, ang):
        _, j = tree.query(p); return np.cos(np.angle(np.exp(1j * (ang - np.angle(Nn[j])))))
    for q, e in enumerate(ev):
        n = e["n"]; ch = dict(e["chain"]); c = ch[n]
        for m in range(n + 1, n + 5):
            h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
            if not len(hp): break
            dd = np.hypot(*(hp - c).T); kk = np.argmin(dd)
            if dd[kk] > 2.5: break
            c = hp[kk]; ch[m] = c
        for lag in LAGS:
            m = n + lag
            if m not in ch: continue
            f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - ch[m]).T))); samp("esc", q, lag, m, ch[m], ctv(ch[m], f[i, 3]))
            if lag % 5 == 0:
                dd = tree.query(f[:, :2])[0]
                res = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (dd < 4)) if np.hypot(*(f[jj, :2] - ch[m])) > 3 and abs(abs(f[jj, 0]) - abs(ch[m][0])) / A < 0.15]
                if res: jj = rng.choice(res); samp("ctl", q, lag, m, f[jj, :2], ctv(f[jj, :2], f[jj, 3]))
    pts = np.c_[XX[inside], YY[inside]]
    for m in range(len(fr) // 2, len(fr), 9):
        f = fr[m]
        if not len(f): continue
        for fd in f: samp("p" if fd[2] > 0 else "m", -1, 0, m, fd[:2], ctv(fd[:2], fd[3]) if fd[2] > 0 else np.nan)
        for pr in pts[rng.choice(len(pts), 25)]: samp("rnd", -1, 0, m, pr + rng.uniform(-0.5, 0.5, 2))
    return out
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    if os.environ.get("REPLOT"):   # redraw from the saved samples (REPLOT=1 SIGMA=2.5 qcrit.py TAG)
        a = np.load(os.path.join(H, f"qcrit_{tag}.npy"), allow_pickle=True)
    else:
      with Pool(min(12, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
      a = np.array(sum(res, []), dtype=object); np.save(os.path.join(H, f"qcrit_{tag}.npy"), a, allow_pickle=True)
      np.save(os.path.join(H, f"qcrit_{tag}_runs.npy"), np.array([(r, i) for r, o in zip(runs, res) for i in range(len(o))], dtype=object), allow_pickle=True)   # run of each row (for the region split)
    kind = a[:, 0].astype(str); lag = a[:, 2].astype(int); d, ct, qn, qd, dq = (a[:, i].astype(float) for i in range(3, 8))
    m_ = lambda x: np.nanmean(x); md = lambda x: np.nanmedian(x); se = lambda x: np.nanstd(x) / np.sqrt(max(np.isfinite(x).sum(), 1))
    L = [f"Q criterion, smoothing sigma {SIG} cell(s). Qn = (|O|^2-|E|^2)/(|O|^2+|E|^2) at the point; d_Q = distance to the nearest Q = 0 line; P(d_Q<1) = fraction within 1 cell",
         "ESCAPERS: lag(t.u.) depth cos_tau  Qn  d_Q(median)  P(d_Q<1)  |  residents: Qn  d_Q  P(d_Q<1)"]
    for lg in LAGS:
        e = (kind == "esc") & (lag == lg); c = (kind == "ctl") & (lag == lg)
        if not e.any(): continue
        s = f"{2*lg:+4d}  {m_(d[e]):.2f}  {m_(ct[e]):+.2f}   {m_(qn[e]):+.3f}±{se(qn[e]):.3f}  {md(dq[e]):.2f}  {np.mean(dq[e] < 1):.2f}"
        if c.any(): s += f"  |  {m_(qn[c]):+.3f}±{se(qn[c]):.3f}  {md(dq[c]):.2f}  {np.mean(dq[c] < 1):.2f}"
        L.append(s)
    DB = [0, 1, 2, 3, 4, 5, 6, 8, 10, 99]; prof = {}
    L.append("\nBY DEPTH: depth | +1/2: n Qn d_Q(med) P(d_Q<1) | -1/2: n Qn d_Q P | random: n Qn d_Q P")
    for i in range(len(DB) - 1):
        row = f"  {DB[i]}-{DB[i+1]}"; dc = 0.5 * (DB[i] + min(DB[i + 1], 14))
        for kd in ("p", "m", "rnd"):
            s_ = (kind == kd) & (d >= DB[i]) & (d < DB[i + 1])
            if s_.sum() >= 20:
                row += f" | {s_.sum():5d} {m_(qn[s_]):+.3f} {md(dq[s_]):.2f} {np.mean(dq[s_] < 1):.2f}"; prof.setdefault(kd, []).append((dc, m_(qn[s_]), md(dq[s_]), np.mean(dq[s_] < 1)))
            else: row += " | -"
        L.append(row)
    L.append("\nBULK (depth > 6), P(d_Q) histogram, ratio to random:  d_Q bin: +1/2  -1/2  random")
    hb = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 10]
    for kd in ("p", "m", "rnd"):
        s_ = (kind == kd) & (d > 6); prof.setdefault("h" + kd, np.histogram(dq[s_], hb, density=True)[0])
    for i in range(len(hb) - 1): L.append(f"  {hb[i]}-{hb[i+1]}: {prof['hp'][i]:.3f}  {prof['hm'][i]:.3f}  {prof['hrnd'][i]:.3f}   (+1/2 / random {prof['hp'][i] / max(prof['hrnd'][i], 1e-9):.2f}, -1/2 / random {prof['hm'][i] / max(prof['hrnd'][i], 1e-9):.2f})")
    print("\n".join(L)); open(os.path.join(H, f"qcrit_{tag}.txt"), "w").write("\n".join(L) + "\n")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4)); tl = 2 * np.array(LAGS)
    # (a) per-event fractions, not the mean Qn: Qn is broadly spread around 0 (median |Qn| 0.3-0.4), so its mean cancels and misleads (user, 2026-10-13).
    #     purple: fraction within 1 cell of a Q = 0 line (the viscometric constraint); orange: fraction with Qn > 0 (rotation side). Solid escaper, dashed resident; bands = binomial s.e.
    cl = sorted(set(lag[kind == "ctl"])); LG = [l for l in LAGS if ((kind == "esc") & (lag == l)).any()]; tl = 2 * np.array(LG)
    def frac(sel, f):
        x = f(sel); n = sel.sum(); p = x.sum() / max(n, 1); return p, np.sqrt(p * (1 - p) / max(n, 1))
    for nm, fn_, col in (("within 1 cell of a Q = 0 line", lambda s_: dq[s_] < 1, "#7a3b8f"), ("Qn > 0 (rotation side)", lambda s_: qn[s_] > 0, "#e08a1e")):
        E = np.array([frac((kind == "esc") & (lag == l), fn_) for l in LG]); C = np.array([frac((kind == "ctl") & (lag == l), fn_) for l in cl])
        ax[0].plot(tl, E[:, 0], "-o", ms=3, color=col, label=f"escaper: {nm}"); ax[0].fill_between(tl, E[:, 0] - E[:, 1], E[:, 0] + E[:, 1], color=col, alpha=0.2, lw=0)
        ax[0].errorbar(2 * np.array(cl), C[:, 0], yerr=C[:, 1], fmt="--s", ms=4, color=col, alpha=0.7, capsize=2, label=f"resident: {nm}")
    ax[0].axhline(np.mean(dq[kind == "rnd"] < 1), color="#7a3b8f", lw=0.8, ls=":", label="random point: within 1 cell")
    ax[0].axhline(0.5, color="0.7", lw=0.6); ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel("time relative to departure (t.u.)"); ax[0].set_ylabel("fraction of samples"); ax[0].axvline(0, color="0.6", lw=0.8, ls=":"); ax[0].set_title("(a) escaper vs resident", fontsize=10)
    for kd, nm, col in (("p", "+1/2", "#c2452d"), ("m", "−1/2", "#1f5fa8"), ("rnd", "random point", "0.5")):
        P_ = np.array(prof[kd]); ax[1].plot(P_[:, 0], P_[:, 1], "-o", ms=3, color=col, label=nm); ax[2].plot(P_[:, 0], P_[:, 3], "-o", ms=3, color=col, label=nm)
    ax[1].set_xlabel("depth (cells from the wall)"); ax[1].set_ylabel("mean Qn at the point"); ax[1].set_title("(b) flow type at defects by depth", fontsize=10)
    ax[2].set_xlabel("depth (cells from the wall)"); ax[2].set_ylabel("P(within 1 cell of a Q = 0 line)"); ax[2].set_title("(c) on viscometric lines?", fontsize=10)
    for a_ in ax[1:]: a_.axhline(0, color="0.8", lw=0.6)
    ax[0].legend(fontsize=6.5, frameon=False, loc="lower left", ncol=1)
    for a_ in ax[1:]: a_.legend(fontsize=7, frameon=False)
    fig.suptitle(f"Q criterion (Head et al. 2024) around defects, smoothing σ = {SIG:g} cell", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_qcrit_{tag}.png"), dpi=200)
