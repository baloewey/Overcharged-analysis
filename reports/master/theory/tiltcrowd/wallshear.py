#!/usr/bin/env python3
"""User's idea (2026-10-15): the wall layer's self-flow opposes the circulation; a denser layer weakens the wall flow more -> less tilt.
Refinement tested here: the tilt law needs simple shear (strain = rotation); the Leslie angle of a flow with strain rate e and rotation rate w/2 obeys
cos 2 theta_L = (w/2) / (lambda e), so tau_pred = arccos(1 / (lambda r_e)) with r_e = |e_rtheta| / |w/2| in the wall band. A counter-flow from the layer can
change r_e (not just the shear magnitude, which the law ignores).
Shape-sweep circles (runs/shapes_server/circle, 45 runs: R_eq 1/1.5/2 l_a, k up to 16, 3 seeds; flowfield every ~? t.u.), second half.
Per flow frame: v_theta(r) azimuthal average (1-cell bins), signed by the frame circulation; e_rtheta = (r/2) d(v_theta/r)/dr, w/2 = (1/2r) d(r v_theta)/dr;
band = 1-4 cells from the wall. Per run: <v_theta>, wall-band strain, rotation, r_e, tau_pred (lambda 2); measured tilt of the coherent layer (C > 0.9) at the
same times; spacing 2 pi R/(2k-2).  -> wallshear.txt, fig_wallshear.png"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); LAM = 2.0
sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/energetics")
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; r = np.hypot(X, Y)
    vth_all = (X * rows[:, 5] - Y * rows[:, 4]) / np.maximum(r, 1e-9)
    edges = np.arange(0, R + 0.01, 1.0); rc = 0.5 * (edges[1:] + edges[:-1])
    ts = ft[len(ft) // 2:]; prof = []
    for tt in ts:
        ii = idx[tt]; m = r[ii] < R; h, _ = np.histogram(r[ii][m], edges, weights=vth_all[ii][m]); c, _ = np.histogram(r[ii][m], edges)
        p = h / np.maximum(c, 1); sgn = np.sign(np.sum(p * rc)); prof.append(p * (sgn if sgn else 1))
    P = np.array(prof); vm = P.mean(0)
    e = 0.5 * rc * np.gradient(vm / rc, rc); w2 = 0.5 * np.gradient(rc * vm, rc) / rc
    band = (rc > R - 4) & (rc < R - 1); bulk = (rc > 0.3 * R) & (rc < R - 6)
    re_b = np.abs(e[band]).mean() / max(np.abs(w2[band]).mean(), 1e-12); re_bulk = np.abs(e[bulk]).mean() / max(np.abs(w2[bulk]).mean(), 1e-12)
    # per-frame r_e in the band (fluctuations)
    t, fr = read_with_angle(run, Q, 0.5); tilts = []
    for f in fr:
        if not len(f): continue
        rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
        if w.sum() < 3: continue
        z = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))
        if abs(z) > 0.9: tilts.append(abs(np.angle(z)))
    tau_pred = np.degrees(np.arccos(min(1.0, 1 / (LAM * re_b)))) if re_b > 0 else np.nan
    return dict(run=os.path.relpath(run, ROOT), k=k, R=R, a=2 * np.pi * R / (2 * k - 2), vth=np.abs(vm[(rc > 0.3 * R) & (rc < R - 2)]).mean(),
                e_band=np.abs(e[band]).mean(), w_band=np.abs(w2[band]).mean(), re_band=re_b, re_bulk=re_bulk, tau_pred=tau_pred,
                tilt=np.degrees(np.median(tilts)) if len(tilts) > 20 else np.nan, prof=vm, rc=rc)


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/shapes_server/circle/circle/Req*/k*/s0*/DONE"))
    if os.environ.get("MAIN10"):   # 10-seed main circle sweep window (k 2-10, R 0.8-2.2 l_a) with fetched flow fields
        runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/flowfield.dat.gz"))
    TAG = os.environ.get("TAG", "")
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs)
    res.sort(key=lambda x: x["a"])
    L = ["WALL-BAND SHEAR vs SPACING (shape-sweep circles, 3 seeds, second half; band 1-4 cells from the wall; lambda 2)",
         "run | k | R | spacing a | |v_theta| | band strain e | band rotation w/2 | r_e band | r_e bulk | tau_pred = arccos(1/(lambda r_e)) | measured tilt (C > 0.9)"]
    for x in res:
        L.append(f"{'/'.join(x['run'].split('/')[-3:]):20s} | {x['k']:2d} | {x['R']:5.1f} | {x['a']:5.1f} | {x['vth']:.3f} | {x['e_band']:.4f} | {x['w_band']:.4f} | "
                 f"{x['re_band']:.2f} | {x['re_bulk']:.2f} | {x['tau_pred']:5.1f} | {x['tilt']:5.1f}")
    a = np.array([x["a"] for x in res]); T = np.array([x["tilt"] for x in res]); P = np.array([x["tau_pred"] for x in res]); RE = np.array([x["re_band"] for x in res])
    ok = np.isfinite(T) & np.isfinite(P)
    L.append(f"corr(measured tilt, spacing) {np.corrcoef(T[ok], np.log(a[ok]))[0,1]:+.2f}; corr(tilt, r_e band) {np.corrcoef(T[ok], RE[ok])[0,1]:+.2f}; corr(tilt, tau_pred) {np.corrcoef(T[ok], P[ok])[0,1]:+.2f}; "
             f"rms(tilt - tau_pred) {np.sqrt(np.mean((T[ok]-P[ok])**2)):.1f} deg")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    sc = ax[0].scatter(a, RE, c=T, cmap="magma", vmin=40, vmax=70, s=18); ax[0].set_xscale("log"); ax[0].set_xlabel("wall-defect spacing (cells)"); ax[0].set_ylabel("strain / rotation in the wall band $r_e$")
    ax[0].axhline(1, color="0.7", ls=":"); fig.colorbar(sc, ax=ax[0], label="measured tilt (deg)"); ax[0].set_title("(a) Wall-band flow vs spacing", fontsize=10)
    ax[1].scatter(P, T, c=np.log(a), cmap="viridis", s=18); lim = [35, 80]; ax[1].plot(lim, lim, "k--", lw=0.8); ax[1].set_xlabel(r"predicted tilt arccos$(1/\lambda r_e)$ (deg)"); ax[1].set_ylabel("measured tilt (deg)")
    ax[1].set_title("(b) Tilt law with the measured strain/rotation ratio", fontsize=10)
    for x in res[::max(1, len(res) // 8)]:
        ax[2].plot(x["R"] - x["rc"], x["prof"] / max(np.abs(x["prof"]).max(), 1e-9), lw=1, label=f"a = {x['a']:.0f}")
    ax[2].set_xlim(0, 12); ax[2].set_xlabel("distance from the wall (cells)"); ax[2].set_ylabel(r"$v_\theta$ / max (signed by circulation)"); ax[2].legend(fontsize=7, frameon=False)
    ax[2].set_title("(c) Near-wall azimuthal flow", fontsize=10)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_wallshear{TAG}.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, f"wallshear{TAG}.txt"), "w").write("\n".join(L) + "\n")
