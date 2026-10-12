#!/usr/bin/env python3
"""Energetics of wall-defect ejection from measured MPCD configurations (user, 2026-10-11): 'the escape is local; what sets the stage may be collective'.
Area-2500 cavities, k = +-6 (runs/fineflow, actlayout, forced_local, swirl_ellipse*). Depleting escapes as in kickout/depleting.py (resident >= 40 t.u., layer <= 2k-2, no swap,
vacancy persists). For each event at frames n (departure), n-10 (20 t.u. before), n-20 (40 t.u. before), and for a resident control in the same frame (same region):
  ELASTIC: all detected defects at their measured positions (frame skipped if their charges do not sum to k); the probe defect moved inward along the wall normal by
  delta = 0, 0.5, 1, 2, 3, 4, 6 cells, others fixed; anisotropic Frank energy (K1 42, K3 23.1) with eMode-2 anchoring, natural relaxed cores (ellipse_binding/binding.py Ellipse).
  -> force F0 = -(E(0.5) - E(0))/0.5 (> 0: pushed inward = out of the layer), barrier = max_delta E(delta) - E(0) along the path, and E(6) - E(0).
  HYDRO: measured flow at the probe (nearest flow frame, cells within 2): inward and tangential components.
Usage: energetics.py RUNDIR [...] -> energetics_<tag>.txt and energetics_<tag>.npz"""
import gzip, json, os, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity"); TH = "/home/baloewe/Physics/Overcharged/reports/master/theory"
sys.path.insert(0, TH + "/ellipse_binding"); os.environ.setdefault("GRID_H", "0.5")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
import binding as bd  # noqa
DELTAS = (0, 0.5, 1, 2, 3, 4, 6); H = os.path.dirname(os.path.abspath(__file__))
def events(run, max_events=12):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; N0 = 2 * abs(k) - 2 if k > 0 else 2 * abs(k) + 2; sg = 1 if k > 0 else -1
    W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max()
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t); D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
    nw = np.array([((np.sign(f[:, 2]) == sg) & (D[n] < 4)).sum() if len(f) else 0 for n, f in enumerate(fr)], float); nws = np.convolve(nw, np.ones(11) / 11, mode="same")
    ev = []; rng = np.random.default_rng(7)
    for n in range(21, len(fr) - 26):
        f, g = fr[n], fr[n + 1]
        if not len(f) or not len(g): continue
        gp = g[np.sign(g[:, 2]) == sg, :2]
        for i in np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] < 4)):
            p = f[i, :2]; dd = np.hypot(*(gp - p).T) if len(gp) else np.array([99.0])
            if dd.min() > 2.5 or tree.query(gp[np.argmin(dd)])[0] < 4: continue
            cur = gp[np.argmin(dd)]; ok = True
            for m in range(n + 2, n + 11):
                h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
                if not len(hp): ok = False; break
                d2 = np.hypot(*(hp - cur).T); kk = np.argmin(d2)
                if d2[kk] > 2.5: ok = False; break
                cur = hp[kk]
            if not ok or tree.query(cur)[0] <= 6: continue
            chain = {n: p}; c = p; res = True
            for m in range(n - 1, n - 21, -1):
                h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
                if not len(hp): res = False; break
                d2 = np.hypot(*(hp - c).T); kk = np.argmin(d2)
                if d2[kk] > 2.5 or tree.query(hp[kk])[0] > 4: res = False; break
                c = hp[kk]; chain[m] = c
            swap = any(len(fr[m]) and ((np.sign(fr[m][:, 2]) == sg) & (D[m] < 4) & (np.hypot(*(fr[m][:, :2] - p).T) < 4) & (np.hypot(*(fr[m][:, :2] - cur).T) > 4)).any() for m in range(n + 1, n + 11))
            if not (nws[n - 5] <= N0 + 0.25 and res and not swap and nws[n + 5:n + 26].max() < nws[n - 5] - 0.25): continue
            pw = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] < 4)) if jj != i and abs(abs(f[jj, 0]) - abs(p[0])) / A < 0.15]
            ev.append(dict(n=n, chain=chain, ctrl=(f[rng.choice(pw), :2] if pw else None)))
    if len(ev) > max_events: ev = [ev[q] for q in np.linspace(0, len(ev) - 1, max_events).astype(int)]
    return J, bc, k, t, fr, ev
def landscape(E, frame, probe, k):
    zs = frame[:, 0] + 1j * frame[:, 1]; qs = np.sign(frame[:, 2]) * 0.5
    if abs(qs.sum() - k) > 1e-6: return None
    ip = int(np.argmin(np.abs(zs - (probe[0] + 1j * probe[1]))))
    nrm = E.Nn[np.argmin(np.abs(E.W - zs[ip]))]; order = [ip] + [j for j in range(len(zs)) if j != ip]
    Es = []
    for dl in DELTAS:
        z = zs[order].copy(); z[0] = zs[ip] - dl * nrm; e, _ = E.energy(list(z), list(qs[order]), [None] * len(z)); Es.append(e)
    Es = np.array(Es) - Es[0]; return dict(F0=-Es[1] / 0.5, barrier=Es.max(), E6=Es[-1], prof=Es)
def flow_at(run, Q, tt, p, cache={}):
    key = run
    if key not in cache:
        fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"; op = gzip.open if fn.endswith(".gz") else open
        rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit())); cache[key] = (rows, np.unique(rows[:, 0]))
    rows, ft = cache[key]; t0 = ft[np.argmin(np.abs(ft - tt))]; r = rows[rows[:, 0] == t0]
    x = r[:, 1] + 0.5 - Q[0]; y = r[:, 2] + 0.5 - Q[1]; m = np.hypot(x - p[0], y - p[1]) < 2.0
    return (r[m, 4].mean(), r[m, 5].mean()) if m.any() else (np.nan, np.nan)
def one(run):
    J, bc, k, t, fr, ev = events(run); eps = bc["aInv"][0] if "aInv" in bc else 1.0; E = bd.Ellipse(float(eps)); Q = bc["Q"]; out = []
    for e in ev:
        for lag in (0, 10, 20):
            m = e["n"] - lag
            if m not in e["chain"]: continue
            for kind, probe in (("escaper", e["chain"][m]), ("control", e["ctrl"] if lag == 0 else None)):
                if probe is None: continue
                L = landscape(E, fr[m], probe, k)
                if L is None: continue
                vx, vy = flow_at(run, Q, t[m], probe); zp = probe[0] + 1j * probe[1]; nrm = E.Nn[np.argmin(np.abs(E.W - zp))]
                vin = -(vx * nrm.real + vy * nrm.imag); vtan = -vx * nrm.imag + vy * nrm.real
                out.append((kind, lag, L["F0"], L["barrier"], L["E6"], vin, abs(vtan)))
    return run, out
if __name__ == "__main__":
    runs = [os.path.abspath(r) for r in sys.argv[1:]]; tag = os.path.basename(os.path.dirname(runs[0]))
    with Pool(min(6, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    L = []; say = lambda x: (L.append(x), print(x, flush=True)); allrows = []
    for run, out in res:
        allrows += [(os.path.basename(run),) + o for o in out]
    for kind, lag in (("escaper", 20), ("escaper", 10), ("escaper", 0), ("control", 0)):
        R = [r for r in allrows if r[1] == kind and r[2] == lag]
        if not R: continue
        a = np.array([r[3:] for r in R], float)
        say(f"{kind:8s} {2 * lag:2d} t.u. before (n {len(R)}): elastic force inward {np.nanmean(a[:, 0]):+7.2f} kBT/cell (median {np.nanmedian(a[:, 0]):+.2f}); barrier {np.nanmean(a[:, 1]):6.1f} "
            f"(median {np.nanmedian(a[:, 1]):.1f}); E(6)-E(0) {np.nanmean(a[:, 2]):+6.1f} kBT; flow inward {np.nanmean(a[:, 3]):+.3f}, |tangential| {np.nanmean(a[:, 4]):.3f}")
    open(os.path.join(H, f"energetics_{tag}.txt"), "w").write("\n".join(L) + "\n")
    np.save(os.path.join(H, f"energetics_{tag}_rows.npy"), np.array(allrows, dtype=object), allow_pickle=True)
