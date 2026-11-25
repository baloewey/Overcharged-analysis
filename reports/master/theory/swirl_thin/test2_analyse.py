#!/usr/bin/env python3
"""Test 2 analysis (user, 2026-11-25): passive layers with the imposed shear confined to the wall (runs/swirl_thin: thin = annulus R-5..R-3, medium = R-9..R-7,
torque-matched) vs the disk forcing (runs/swirl_crowd A03_k*, swirl 0.03 within 0.6R; disk_k04 in swirl_thin). Per run, second half:
  flow profile <v_theta>(y) (flow every 10 t.u.) -> strain rate |e|(y) = |r d(v/r)/dr|/2 and rotation; fraction of the strain within the outer 3 / 7 cells (check)
  wall +1/2 count (2k vs 2k-2: ejection), bulk charge, layer tilt over coherent frames (C > 0.9) and coherent fraction (chi1_compare.one)
  band profile beyond the near field (static-subtracted; fig_band_profiles_k.run_one) at 2.5 / 4.5 / 6.5 cells
-> test2.txt, test2.npy"""
import glob, gzip, json, os, sys, numpy as np
from multiprocessing import Pool
TH = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path[:0] = [TH + "/chi1"]
import chi1_compare as CC, fig_band_profiles_k as BP  # noqa
H = os.path.dirname(os.path.abspath(__file__)); RUNS = "/home/baloewe/Physics/Overcharged/runs"


def flow_profile(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; Q = np.array(bc["Q"][:2], float)
    fn = run + "/flowfield.dat.gz"; op = gzip.open if os.path.exists(fn) else open; fn = fn if os.path.exists(fn) else run + "/flowfield.dat"
    with op(fn, "rt") as fh: A = np.loadtxt((l for l in fh if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    ts = np.unique(A[:, 0]); A = A[A[:, 0] >= ts[len(ts) // 2]]
    x = A[:, 1] + 0.5 - Q[0]; y = A[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); vt = (x * A[:, 5] - y * A[:, 4]) / np.maximum(r, 1e-9)
    edges = np.r_[np.arange(0, R, 1.0), R]; b = np.digitize(r, edges) - 1; ok = (r < R) & (b >= 0) & (b < len(edges) - 1)
    v = np.bincount(b[ok], vt[ok], len(edges) - 1) / np.maximum(np.bincount(b[ok], None, len(edges) - 1), 1); rc = 0.5 * (edges[1:] + edges[:-1])
    e = np.abs(0.5 * rc * np.gradient(v / rc, rc)); dep = R - rc; m = dep > 0.5
    return dict(e_out3=e[m & (dep < 3.5)].sum() / e[m].sum(), e_out7=e[m & (dep < 7.5)].sum() / e[m].sum(), vmax=v.max(), s_wall=np.interp(1.0, dep[::-1], (e * 2)[::-1]),
                prof=(dep[::-1], v[::-1]))


def one(run):
    try:
        a = CC.one(run); b = BP.run_one(run); f = flow_profile(run)
        band = {y: (np.interp(y, BP.YB + 0.5, b["prof"]) if b and "err" not in b else np.nan) for y in (2.5, 4.5, 6.5)}
        name = os.path.basename(run); cond = "disk" if name.startswith(("A03", "disk")) else name.split("_")[0]
        return dict(run=name, cond=cond, k=a["k"], nw=a["nw"], qb=a["qb"], tilt=a["tilt"], fco=a["fco"], edge=b["edge"] if b and "err" not in b else np.nan, band=band, **{k_: f[k_] for k_ in ("e_out3", "e_out7", "vmax", "s_wall")})
    except Exception as ex:
        return dict(run=run, err=repr(ex)[:150])


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(RUNS + "/swirl_thin/*/DONE")) + sorted(d for d in glob.glob(RUNS + "/swirl_crowd/A03_k*_s*") if os.path.exists(d + "/DONE") or os.path.exists(d + "/directorfield.dat.gz"))
    with Pool(12) as p: res = p.map(one, runs, chunksize=1)
    ok = [r for r in res if "err" not in r]
    L = ["TEST 2: shear confined to the wall (thin: outer 3 cells; medium: outer 7) vs disk forcing (A03), passive, R 31.83, lambda 2, torque-matched. Means over seeds.",
         "cond | k | runs | strain in outer 3 / 7 cells | wall shear | wall +1/2 (2k-2) | bulk | tilt (coherent frac) | 2 x band at 2.5 / 4.5 / 6.5 cells (static-subtracted)"]
    for e in [r for r in res if "err" in r]: L.append(f"  ERROR {e['run'][-40:]}: {e['err']}")
    for cond in ("disk", "medium", "thin"):
        for k in (4, 6, 8, 10, 13, 16):
            g = [r for r in ok if r["cond"] == cond and r["k"] == k]
            if not g: continue
            f = lambda key: np.nanmean([r[key] for r in g]); bd = lambda y: 2 * np.nanmean([r["band"][y] for r in g])
            L.append(f"{cond:6s} | {k:2d} | {len(g)} | {f('e_out3'):.2f} / {f('e_out7'):.2f} | {f('s_wall'):.4f} | {f('nw'):5.2f} ({2 * k - 2}) | {f('qb'):+.2f} | {f('tilt'):5.1f} ({f('fco'):.2f}) | "
                     f"{bd(2.5):5.1f} / {bd(4.5):5.1f} / {bd(6.5):5.1f}")
    np.save(os.path.join(H, "test2.npy"), ok, allow_pickle=True); open(os.path.join(H, "test2.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
