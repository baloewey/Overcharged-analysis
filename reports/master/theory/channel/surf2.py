#!/usr/bin/env python3
"""Fast pockets and surfing +1/2, part 2 (user, 2026-10-19).
(A) BEND-WALL TEST at director frames (every 100 t.u., t >= 3000): active force density |div Q| (Q = S(nn - I/2), centred differences) and the bend |(n.grad)n| and
    splay |div n| magnitudes (director signs aligned to each cell's own n before differencing); at the fastest 2 % of inner-channel cells vs all inner cells
    (ratio of means), Spearman correlation of speed with |div Q| (smoothed 2 cells), and the fraction of fast cells within 2 cells of the top-10 % |div Q| cells.
(B) SURFING EVENTS (window 3000-4000 t.u., defects and flow every 2 t.u.): bulk +1/2 tracks (nearest-neighbour links < 2.5 cells, periodic in x); 'surfing'
    segments = displacement >= 25 cells within 40 t.u. (mean speed >= 0.6); per segment: duration, length, <cos(motion, head)>, <cos(motion, local flow)>,
    local flow speed / channel rms, and the head-flow alignment <cos(head, local flow)>; compared with all bulk +1/2 steps. -> surf2.txt"""
import os, sys, json, gzip, numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
def frames(fn, tmin, tmax, ncol):
    F = {}
    with gzip.open(fn, "rt") as fh:
        for l in fh:
            p = l.split()
            if len(p) != ncol or not p[0][0].isdigit(): continue
            tt = float(p[0])
            if tt < tmin: continue
            if tt > tmax: break
            F.setdefault(tt, []).append([float(v) for v in p[1:]])
    return {k: np.array(v) for k, v in F.items()}
def grid(a, L, Hd, cols):
    out = [np.zeros((L, Hd)) for _ in cols]; ix, iy = a[:, 0].astype(int), a[:, 1].astype(int)
    for o, c in zip(out, cols): o[ix, iy] = a[:, c]
    return out
def partA(run, L, Hd, y0, y1):
    D = frames(run + "/directorfield.dat.gz", 3000, 6000, 8); Fl = frames(run + "/flowfield.dat.gz", 3000, 6000, 7); yc = np.arange(Hd) + 0.5; inner = (yc > y0 + 3) & (yc < y1 - 3)
    rq, rb, rs, rho, near = [], [], [], [], []
    for t in sorted(D):
        if t not in Fl: continue
        nx, ny, S = grid(D[t], L, Hd, (3, 4, 6)); U, V = grid(Fl[t], L, Hd, (3, 4))
        Qxx, Qxy = S * (nx * nx - 0.5), S * nx * ny
        dQxx_x = (np.roll(Qxx, -1, 0) - np.roll(Qxx, 1, 0)) / 2; dQxy_y = np.gradient(Qxy, axis=1); dQxy_x = (np.roll(Qxy, -1, 0) - np.roll(Qxy, 1, 0)) / 2; dQxx_y = np.gradient(Qxx, axis=1)
        fx, fy = dQxx_x + dQxy_y, dQxy_x - dQxx_y; fq = gaussian_filter(np.hypot(fx, fy), 1.0, mode="wrap")
        def aligned(dx, dy):     # neighbour director with sign matched to the local one
            mx, my = np.roll(np.roll(nx, -dx, 0), -dy, 1), np.roll(np.roll(ny, -dx, 0), -dy, 1); s = np.sign(mx * nx + my * ny); s[s == 0] = 1; return mx * s, my * s
        exp_, eyp = aligned(1, 0); exm, eym = aligned(-1, 0); ezp, ezq = aligned(0, 1); ewm, ewn = aligned(0, -1)
        dnx_dx, dny_dx = (exp_ - exm) / 2, (eyp - eym) / 2; dnx_dy, dny_dy = (ezp - ewm) / 2, (ezq - ewn) / 2
        splay = np.abs(dnx_dx + dny_dy); bx = nx * dnx_dx + ny * dnx_dy; by = nx * dny_dx + ny * dny_dy; bend = np.hypot(bx, by)
        sp = np.hypot(U, V); m = np.zeros_like(sp, bool); m[:, inner] = True; thr = np.percentile(sp[m], 98); fast = m & (sp >= thr)
        rq.append(fq[fast].mean() / fq[m].mean()); rb.append(bend[fast].mean() / bend[m].mean()); rs.append(splay[fast].mean() / splay[m].mean())
        rho.append(spearmanr(sp[m], fq[m])[0]); top = m & (fq >= np.percentile(fq[m], 90)); topd = gaussian_filter(top.astype(float), 0.01) > 0
        from scipy.ndimage import binary_dilation; near.append(binary_dilation(top, iterations=2)[fast].mean())
    return np.mean(rq), np.mean(rb), np.mean(rs), np.mean(rho), np.mean(near), len(rq)
def partB(run, L, Hd, y0, y1, T0=3000.0, TL=1000.0):
    Fl = frames(run + "/flowfield.dat.gz", T0, T0 + TL, 7); t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); sel = np.flatnonzero((t >= T0) & (t <= T0 + TL))
    rms = np.sqrt(np.mean([np.mean(a[:, 3] ** 2 + a[:, 4] ** 2) for a in Fl.values()]))
    tracks, prev = [], []
    for i in sel:
        f = fr[i]; p = f[(f[:, 2] > 0) & (f[:, 1] - y0 >= 4) & (y1 - f[:, 1] >= 4)] if len(f) else np.zeros((0, 4)); new = []
        used = set()
        for x, y, q, a in p:
            best, bd = None, 2.5
            for tr in prev:
                if id(tr) in used: continue
                lx, ly = tr[-1][1], tr[-1][2]; d = np.hypot((x - lx + L / 2) % L - L / 2, y - ly)
                if d < bd: best, bd = tr, d
            if best is not None: used.add(id(best)); best.append((t[i], x, y, a)); new.append(best)
            else: tr = [(t[i], x, y, a)]; tracks.append(tr); new.append(tr)
        prev = new
    def flow_at(tt, x, y):
        a = Fl.get(tt)
        if a is None: return np.array([np.nan, np.nan])
        m = (np.abs((a[:, 0] + 0.5 - x + L / 2) % L - L / 2) < 2) & (np.abs(a[:, 1] + 0.5 - y) < 2); return a[m, 3:5].mean(0) if m.any() else np.array([np.nan, np.nan])
    seg, allst = [], []
    for tr in tracks:
        tr = np.array(tr)
        if len(tr) < 3: continue
        xs = np.unwrap(tr[:, 1] * 2 * np.pi / L) * L / (2 * np.pi); ys = tr[:, 2]
        for j in range(len(tr) - 1):
            mv = np.array([xs[j + 1] - xs[j], ys[j + 1] - ys[j]]); h = np.array([np.cos(tr[j, 3] + np.pi), np.sin(tr[j, 3] + np.pi)]); u = flow_at(tr[j, 0], tr[j, 1], tr[j, 2])
            allst.append((mv, h, u))
        n20 = 20                                     # 40 t.u. at 2 t.u. per frame
        j = 0
        while j + n20 < len(tr):
            if np.hypot(xs[j + n20] - xs[j], ys[j + n20] - ys[j]) >= 25:
                k = j + n20
                while k + 1 < len(tr) and np.hypot(xs[k + 1] - xs[k], ys[k + 1] - ys[k]) > 0.4: k += 1
                st = [allst[-(len(tr) - 1) + m] for m in range(j, k)] if len(allst) >= len(tr) - 1 else []
                seg.append(dict(dur=tr[k, 0] - tr[j, 0], length=np.hypot(xs[k] - xs[j], ys[k] - ys[j]), steps=st)); j = k
            j += 1
    cos = lambda a, b: float(np.dot(a, b) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))
    def stats(st):
        st = [s for s in st if np.all(np.isfinite(s[2])) and np.linalg.norm(s[0]) > 0.05]
        if not st: return (np.nan,) * 4
        return (np.mean([cos(m, h) for m, h, u in st]), np.mean([cos(m, u) for m, h, u in st]), np.mean([np.linalg.norm(u) for m, h, u in st]) / rms, np.mean([cos(h, u) for m, h, u in st]))
    allS = stats(allst); segS = stats([s for g in seg for s in g["steps"]])
    return len(tracks), len(seg), np.median([g["dur"] for g in seg]) if seg else np.nan, np.median([g["length"] for g in seg]) if seg else np.nan, allS, segS, rms
def one(r):
    run = f"{R_}/{r}"; J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    A = partA(run, L, Hd, y0, y1); B = partB(run, L, Hd, y0, y1)
    return [f"{r}",
            f"  (A) {A[5]} director frames: at the fastest 2% of cells |div Q| is {A[0]:.2f}x its mean, bend {A[1]:.2f}x, splay {A[2]:.2f}x; Spearman(speed, |div Q|) {A[3]:+.2f}; fast cells within 2 cells of the top-10% |div Q|: {A[4]:.2f}",
            f"  (B) {B[0]} bulk +1/2 tracks; surfing segments (>= 25 cells in 40 t.u.): {B[1]}, median duration {B[2]:.0f} t.u., median length {B[3]:.0f} cells",
            f"      all steps: cos(motion, head) {B[4][0]:+.2f}, cos(motion, local flow) {B[4][1]:+.2f}, local flow / rms {B[4][2]:.2f}, cos(head, local flow) {B[4][3]:+.2f}",
            f"      surfing:   cos(motion, head) {B[5][0]:+.2f}, cos(motion, local flow) {B[5][1]:+.2f}, local flow / rms {B[5][2]:.2f}, cos(head, local flow) {B[5][3]:+.2f}"]
if __name__ == "__main__":
    runs = ["channel_current/W24_sym_s5", "channel_current/W24_mixed_s3", "channel_pressure/sym_g+0.001_s2", "channel_pressure/mixed_g+0.003_s1"]
    with Pool(4) as p: res = p.map(one, runs)
    L = ["FAST POCKETS AND SURFING +1/2, part 2 (surf2.py)"] + [l for r in res for l in r]
    open(os.path.join(H, "surf2.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
