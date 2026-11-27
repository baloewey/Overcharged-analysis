#!/usr/bin/env python3
"""Lock slope with SUB-CELL core positions (user, 2026-11-27): removes the regression dilution from the trackers' grid-rounded positions (+-0.5 cell,
sd 0.29). For each director frame (second half) whose matching defect frame has exactly 2k-2 wall +1/2: build q1 = S cos 2theta, q2 = S sin 2theta on the
cell centres, with the cells within 1.5 of the wall outside filled with the anchoring (theta = k phi, S = 1, the imposed BC); for each tracked wall +1/2 find the nearest plaquette (2x2 cell centres) within 2.5 cells with winding +1/2 and solve
the bilinear interpolant q1 = q2 = 0 inside it (Newton; fallback: plaquette centre). Then the fluctuation regression of lock_slope_fluct.py (dtau_j vs
dphi_j about the frame's regular ring) with (a) the tracker's positions and (b) the refined positions, same frames, same tilts (tracker angle).
Sets: active circles R 31.83 (main, director every 100 t.u.), Test 2 thin (every 20 t.u.), disk A03 (every 100 t.u.). -> lock_slope_subcell.txt"""
import glob, gzip, json, os, sys, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa
RR = "/home/baloewe/Physics/Overcharged/runs"


def wind(a):                                      # headless winding of director angles around a closed loop (list of angles), in units of 2 pi
    d = np.diff(np.r_[a, a[0]]); d = (d + np.pi / 2) % np.pi - np.pi / 2; return d.sum() / (2 * np.pi)


def refine(q1, q2, th, ok, x0, y0):
    """sub-cell zero of Q near grid-index position (x0, y0) (cell-centre index coordinates); None if not found"""
    best = None
    for i in range(int(np.floor(x0)) - 3, int(np.floor(x0)) + 3):
        for j in range(int(np.floor(y0)) - 3, int(np.floor(y0)) + 3):
            if i < 0 or j < 0 or i + 1 >= ok.shape[0] or j + 1 >= ok.shape[1] or not ok[i:i + 2, j:j + 2].all(): continue
            if abs(wind([th[i, j], th[i + 1, j], th[i + 1, j + 1], th[i, j + 1]]) - 0.5) > 1e-6: continue
            d = np.hypot(i + 0.5 - x0, j + 0.5 - y0)
            if d < 2.5 and (best is None or d < best[0]): best = (d, i, j)
    if best is None: return None
    _, i, j = best; A = np.array([[q1[i, j], q1[i + 1, j], q1[i, j + 1], q1[i + 1, j + 1]], [q2[i, j], q2[i + 1, j], q2[i, j + 1], q2[i + 1, j + 1]]])
    u = v = 0.5
    for _ in range(30):
        f = A @ np.array([(1 - u) * (1 - v), u * (1 - v), (1 - u) * v, u * v]); J = A @ np.array([[-(1 - v), -(1 - u)], [1 - v, -u], [-v, 1 - u], [v, u]])
        try: du = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError: break
        u, v = u + du[0], v + du[1]
        if np.abs(du).max() < 1e-8: break
    if not (-0.05 <= u <= 1.05 and -0.05 <= v <= 1.05): u = v = 0.5
    return i + np.clip(u, 0, 1), j + np.clip(v, 0, 1)


def run_one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = 2 * k - 2
        with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; t, fr = SB.read_with_angle(run, np.zeros(2), 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; out = {"grid": ([], []), "sub": ([], []), "grid_all": ([], []), "sub_all": ([], [])}; nref = ntot = 0
        for tt in td:
            i0 = int(np.argmin(np.abs(t - tt)))
            if abs(t[i0] - tt) > 1e-6: continue
            f = fr[i0]; zs = (f[:, 0] - Q[0]) + 1j * (f[:, 1] - Q[1]); w = (f[:, 2] > 0) & (R - np.abs(zs) < 4)
            if w.sum() != N: continue
            dd = D[D[:, 0] == tt]; th = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th[ix, iy] = np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]
            X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij"); rr_ = np.hypot(X, Y)
            ghost = (rr_ >= R - 0.5) & (rr_ < R + 1.5); th = np.where(ghost, k * np.arctan2(Y, X), th); S = np.where(ghost, 1.0, S)   # wall cells: anchoring theta = k phi (lc.c), fully ordered
            ok = np.isfinite(th) & (rr_ < R + 1.5)
            q1, q2 = S * np.cos(2 * np.nan_to_num(th)), S * np.sin(2 * np.nan_to_num(th))
            pos_g, pos_s = [], []; nr0 = nref
            for xf, yf in zip(f[w, 0], f[w, 1]):
                ntot += 1; rr = refine(q1, q2, th, ok, xf - 0.5, yf - 0.5)          # tracker coordinates -> cell-centre index coordinates
                pos_g.append((xf - Q[0]) + 1j * (yf - Q[1]))
                if rr is None: pos_s.append(pos_g[-1])
                else: nref += 1; pos_s.append((rr[0] + 0.5 - Q[0]) + 1j * (rr[1] + 0.5 - Q[1]))
            tails = f[w, 3]; allref = (nref - nr0) == N
            for key, P in (("grid", np.array(pos_g)), ("sub", np.array(pos_s))) + ((("grid_all", np.array(pos_g)), ("sub_all", np.array(pos_s))) if allref else ()):
                phi = np.angle(P); o = np.argsort(phi); phi = phi[o]; tau = np.angle(np.exp(1j * (tails[o] - phi)))      # tilt relative to the normal at the used position
                alpha = np.angle(np.mean(np.exp(1j * N * phi))) / N; j = np.round((phi - alpha) / (2 * np.pi / N))
                dphi = np.angle(np.exp(1j * (phi - alpha - 2 * np.pi * j / N))); dtau = np.angle(np.exp(1j * (tau - np.angle(np.mean(np.exp(1j * tau))))))
                out[key][0].extend(dphi * np.abs(P).mean()); out[key][1].extend(dtau)
        return dict(k=k, R=R, out=out, frac=nref / max(ntot, 1), run=run)
    except Exception as ex:
        return dict(err=repr(ex)[:160], run=run)


if __name__ == "__main__":
    main = [os.path.dirname(p) for p in sorted(glob.glob(RR + "/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
            if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    sets = [("Active circles", main), ("Driven: shear in outer 3 cells", sorted(glob.glob(RR + "/swirl_thin/thin_k*_s*"))),
            ("Driven: disk A03", sorted(glob.glob(RR + "/swirl_crowd/A03_k*_s*")))]
    allr = [r for _, rs in sets for r in rs]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = dict(zip(allr, p.map(run_one, allr, chunksize=1)))
    L = ["LOCK SLOPE WITH SUB-CELL CORE POSITIONS (zero of the bilinear Q interpolant) vs tracker grid positions; same frames, same tilts",
         "set | k | defect-frames | refined fraction | sd displacement grid / sub (cells) | slope/(2k-2): grid -> sub | corr grid -> sub || same, only frames with ALL wall defects refined"]
    for lab, rs in sets:
        ok = [res[r] for r in rs if res[r] and "err" not in res[r]]
        for e in [res[r] for r in rs if res[r] and "err" in res[r]][:2]: L.append(f"  ERROR {e['run'][-40:]}: {e['err']}")
        for k in sorted({r["k"] for r in ok}):
            g = [r for r in ok if r["k"] == k]; row = {}
            for key in ("grid", "sub", "grid_all", "sub_all"):
                x = np.concatenate([np.array(r["out"][key][0]) for r in g]); y = np.concatenate([np.array(r["out"][key][1]) for r in g])
                if len(x) < 100: row[key] = (len(x), np.nan, np.nan, np.nan); continue
                xr = x / np.mean([r["R"] for r in g]); row[key] = (len(x), x.std(), np.cov(xr, y)[0, 1] / np.var(xr) / (2 * k - 2), np.corrcoef(x, y)[0, 1])
            L.append(f"{lab:32s} | {k:2d} | {row['grid'][0]:6d} | {np.mean([r['frac'] for r in g]):.2f} | {row['grid'][1]:.2f} / {row['sub'][1]:.2f} | "
                     f"{row['grid'][2]:.2f} -> {row['sub'][2]:.2f} | {row['grid'][3]:.2f} -> {row['sub'][3]:.2f} || all-refined frames: n {row['grid_all'][0]}, sd {row['grid_all'][1]:.2f} / {row['sub_all'][1]:.2f}, slope {row['grid_all'][2]:.2f} -> {row['sub_all'][2]:.2f}")
    open(os.path.join(H, "lock_slope_subcell.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
