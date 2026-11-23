#!/usr/bin/env python3
"""Flow-at-core test (pre-registered 2026-10-28 in RESEARCH.md; analysed 2026-11-23). Does crowding weaken the hydrodynamic drive off the wall (a smaller
inward flow at wall cores) or strengthen the elastic hold (the defect moving less than the flow carries it)?
Data: MPCD runs/coreflow (aspect-2 ellipse R_eq 1.5 l_a, act 0.07, k 5 / 14, 5 seeds, flow averaged over 2 t.u. = the defect cadence, 5000 t.u., first 10 %
dropped); continuum runs/continuum_Ktest/cf{1,8}_k*_s* (ellipse_active_K.py --coreflow 1, zeta 1 / 8, 2 seeds, t >= 500; flow per defect row).
Per wall +1/2 (validated tracker track_new, gate 9; wall: distance d < 4), n = inward wall normal at the nearest wall point:
  u_n  flow at the core . n (disk r < 2 around the core, cells inside the cavity; MPCD: mean of the flow frames t and t + 2, i.e. averaged over t +- 2)
  v_n  defect velocity . n, central difference over +-2 t.u. (frames n-1, n+1)
  w_n  v_n - u_n (elastic forces + self-propulsion: the motion relative to the fluid)
  u_t  |flow along the wall|; MPCD only: cos tau = head direction . n (head = +1/2 comet head, from the defect angle)
Sets: all wall frames; residents (d < 4 for the next 20 t.u.); pre-departure (the 10 t.u. before the last wall frame of a track that then goes deeper than 6
within 20 t.u.). Predictions: (a) stronger hold when crowded: u_n similar or larger, w_n more negative; (b) weaker drive: u_n smaller, w_n similar.
-> coreflow.txt, coreflow_rows.npz"""
import glob, gzip, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); TH = os.path.join(H, "..")
sys.path[:0] = [os.path.join(TH, "crowd_protect", "tracker_validation"), "/home/baloewe/Physics/Overcharged/runs/cavity"]
from trackers import track_new  # noqa
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
R = "/home/baloewe/Physics/Overcharged/runs"


class Wall:
    def __init__(self, W):
        self.W = W; self.tree = cKDTree(W); t = np.roll(W, -1, 0) - np.roll(W, 1, 0); t /= np.hypot(*t.T)[:, None]
        n = np.column_stack([-t[:, 1], t[:, 0]]); n *= np.sign(np.sum(n * (-W), 1))[:, None]; self.n, self.t = n, t   # inward normal, tangent
        self.A, self.B = W[:, 0].max(), W[:, 1].max()
    def q(self, xy): d, j = self.tree.query(xy); return d, self.n[j], self.t[j]
    def inside(self, x, y): return (x / self.A) ** 2 + (y / self.B) ** 2 < 1


def mpcd_flow(run, nx, ny):
    with gzip.open(run + "/flowfield.dat.gz", "rb") as f: raw = f.read()
    i = raw.index(b"\n 0.00000e+00\t") + 1; a = np.fromstring(raw[i:].decode(), dtype=np.float32, sep=" ").reshape(-1, 7); del raw
    nf = len(a) // (nx * ny); a = a[:nf * nx * ny]; t = a[::nx * ny, 0].astype(float)
    return t, a[:, 4].reshape(nf, nx, ny), a[:, 5].reshape(nf, nx, ny)


def analyse(fr, t, flow_at, wall, ang=None):
    """fr: frames (x, y, q, ...) centred; flow_at(frame index, x, y) -> (ux, uy); returns rows (set, u_n, v_n, w_n, u_t, cos_tau, d)."""
    rows = []; nf = len(fr)
    for tr in (tr for tr in track_new([f[:, :3] for f in fr], 9.0) if tr["sg"] > 0):
        nn = np.array(tr["n"]); xy = np.array(tr["xy"]); d, nrm, tg = wall.q(xy); pos = dict(zip(nn, range(len(nn))))
        dep = set()
        for i, m in enumerate(nn):                                      # departures: last wall frame before going deeper than 6 within 20 t.u.
            if d[i] < 4 and any(j in pos and d[pos[j]] > 6 for j in range(m + 1, m + 11)):
                j6 = next(j for j in range(m + 1, m + 11) if j in pos and d[pos[j]] > 6)
                if all(j not in pos or d[pos[j]] >= 4 for j in range(m + 1, j6)): dep.add(m)
        pre = {m - l for m in dep for l in range(0, 6)}
        for i, m in enumerate(nn):
            if d[i] >= 4 or (m - 1) not in pos or (m + 1) not in pos or m + 1 >= nf: continue
            v = (xy[pos[m + 1]] - xy[pos[m - 1]]) / (t[m + 1] - t[m - 1]); u = flow_at(m, *xy[i])
            if not np.all(np.isfinite(u)): continue
            res = all(j in pos and d[pos[j]] < 4 for j in range(m + 1, m + 11)); s = 2 if m in pre else (1 if res else 0)
            ct = np.nan
            if ang is not None:
                f = fr[m]; k = np.argmin(np.hypot(f[:, 0] - xy[i, 0], f[:, 1] - xy[i, 1])); a = f[k, 3]; ct = -(np.cos(a) * nrm[i, 0] + np.sin(a) * nrm[i, 1])   # angle column = tail direction (depleting.py convention): head = -angle
            un, vn = u @ nrm[i], v @ nrm[i]; rows.append((s, un, vn, vn - un, abs(u @ tg[i]), ct, d[i]))
    return np.array(rows)


def one_mpcd(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = np.array(bc["Q"][:2], float); W, _, _ = wall_polyline(bc); wall = Wall(W)
    t, fr = read_with_angle(run, Q, 0.1); tf, UX, UY = mpcd_flow(run, *J["domain"])
    X, Y = np.meshgrid(np.arange(J["domain"][0]) + 0.5 - Q[0], np.arange(J["domain"][1]) + 0.5 - Q[1], indexing="ij"); IN = wall.inside(X, Y)
    fi = {round(x, 3): i for i, x in enumerate(tf)}
    def flow_at(m, x, y):
        i = fi.get(round(t[m], 3))
        if i is None or i + 1 >= len(tf): return np.array([np.nan, np.nan])
        c = (np.hypot(X - x, Y - y) < 2.0) & IN
        return np.array([(UX[i][c].mean() + UX[i + 1][c].mean()) / 2, (UY[i][c].mean() + UY[i + 1][c].mean()) / 2]) if c.any() else np.array([np.nan, np.nan])
    return analyse(fr, t, flow_at, wall, ang=True)


def one_cont(base):
    Z = np.load(base + "_defects.npz"); t, n, rows = Z["t"], Z["n"], Z["rows"]; off = np.r_[0, np.cumsum(n)]; keep = t >= 500
    fr = [rows[off[i]:off[i + 1]] for i in range(len(t)) if keep[i]]; t = t[keep]
    Req = 31.830989; th = np.linspace(0, 2 * np.pi, 8000, endpoint=False); wall = Wall(np.column_stack([Req * np.sqrt(2) * np.cos(th), Req / np.sqrt(2) * np.sin(th)]))
    def flow_at(m, x, y):
        f = fr[m]; k = np.argmin(np.hypot(f[:, 0] - x, f[:, 1] - y)); return f[k, 3:5]
    return analyse(fr, t, flow_at, wall)


def summary(lab, rr):
    out = []
    for name, sel in (("all wall", lambda r: r), ("residents", lambda r: r[r[:, 0] == 1]), ("pre-departure", lambda r: r[r[:, 0] == 2])):
        x = sel(rr)
        if not len(x): out.append(f"{lab:14s} {name:13s} n 0"); continue
        se = lambda c: x[:, c].std() / np.sqrt(len(x))
        out.append(f"{lab:14s} {name:13s} n {len(x):6d} | u_n {x[:, 1].mean():+.4f} +- {se(1):.4f} | v_n {x[:, 2].mean():+.4f} +- {se(2):.4f} | "
                   f"w_n {x[:, 3].mean():+.4f} +- {se(3):.4f} | |u_t| {x[:, 4].mean():.4f} | cos tau {np.nanmean(x[:, 5]) if np.isfinite(x[:, 5]).any() else np.nan:+.3f} | d {x[:, 6].mean():.2f}")
    return out


if __name__ == "__main__":
    mp = sorted(glob.glob(R + "/coreflow/ell2_k*_s*")); cb = sorted(f[:-len("_defects.npz")] for f in glob.glob(R + "/continuum_Ktest/cf*_k*_s*_defects.npz"))
    with Pool(3) as p: rm = p.map(one_mpcd, mp)
    with Pool(8) as p: rc = p.map(one_cont, cb)
    L = ["FLOW AT WALL CORES (pre-registered 2026-10-28). Means +- s.e. over core-frames (frames are correlated: s.e. are lower bounds; seed spread below).",
         "u_n inward flow at the core; v_n inward defect velocity; w_n = v_n - u_n; |u_t| flow along the wall; cos tau head . inward normal (MPCD); d depth"]
    S = {}
    for lab, runs, res in (("MPCD", mp, rm), ("cont", cb, rc)):
        for run, r in zip(runs, res):
            b = os.path.basename(run); key = (lab + " " + b.split("_s")[0].replace("ell2_", "")); S.setdefault(key, []).append((b, r))
    for key, v in S.items():
        L += summary(key, np.concatenate([r for _, r in v if len(r)]))
        L.append(f"{'':14s} seed means u_n: " + ", ".join(f"{r[:, 1].mean():+.4f}" for _, r in v if len(r)) + " | w_n: " + ", ".join(f"{r[:, 3].mean():+.4f}" for _, r in v if len(r)))
    np.savez_compressed(os.path.join(H, "coreflow_rows.npz"), **{k.replace(" ", "_"): np.concatenate([r for _, r in v if len(r)]) for k, v in S.items()})
    open(os.path.join(H, "coreflow.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
