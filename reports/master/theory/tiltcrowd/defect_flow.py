#!/usr/bin/env python3
"""Defect-centred wall flow (user's idea, 2026-10-15): do the self-flow lobes of close neighbours change the flow AT a wall +1/2, and with it its tilt?
(The azimuthal mean of the wall flow does not change with spacing, wallshear.txt; the lobes are periodic with the spacing and cancel in that mean.)
Shape-sweep circles with flowfield (runs/shapes_server/circle; 45 runs; flow every 100 t.u. = interval averages; second half).
Per flow frame: wall +1/2 (within 4 cells; positions from the defect frame nearest the flow time). Around each, the flow (cells, bilinear interpolation) is
sampled on a local grid: x along the wall in the circulation direction, y = depth into the cavity; averaged over all wall defects and frames of a run.
At the mean defect depth: velocity gradient by centred differences (+-1 cell) -> strain tensor E (e_xx, e_xy), rotation W = (du_x/dy - du_y/dx)/2
(> 0 for wall shear in the circulation sense), |E| = sqrt(e_xx^2 + e_xy^2), r = |E|/W, phi = atan2(e_xx, e_xy) (0 = simple shear along the wall).
Leslie fixed point of Jeffery's equation in a general 2D flow: tau_pred = arccos(W/(lambda |E|)) - phi  (simple shear: arccos(1/lambda)); reported with
the measured tilt (coherent frames, C > 0.9) and its components. Same at depths 1-6 cells (profile).  -> defect_flow.txt, fig_defect_flow.png"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.ndimage import map_coordinates
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); LAM = 2.0
sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/energetics")
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
XS = np.arange(-6, 6.01, 0.5); YS = np.arange(-1, 8.01, 0.5)      # local grid: along the wall, depth relative to the defect


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    rows, ft, idx = load_flow(run); nx, ny = rows[:, 1].astype(int).max() + 1, rows[:, 2].astype(int).max() + 1
    t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
    acc = np.zeros((2, len(XS), len(YS))); nacc = 0; depths = []; tilts = []
    GX, GY = np.meshgrid(XS, YS, indexing="ij")
    for tt in ft[len(ft) // 2:]:
        ii = idx[tt]; U = np.zeros((nx, ny)); V = np.zeros((nx, ny)); U[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 4]; V[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 5]
        cx = np.arange(nx) + 0.5 - Q[0]; cy = np.arange(ny) + 0.5 - Q[1]; XX, YY = np.meshgrid(cx, cy, indexing="ij"); rr = np.hypot(XX, YY)
        circ = np.sign(np.sum((XX * V - YY * U)[rr < R - 1]))
        if circ == 0: continue
        f = fr[int(np.argmin(np.abs(t - tt)))]
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        if w.sum() < 3: continue
        z = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))
        if abs(z) > 0.9: tilts.append(circ * np.angle(z))
        for x0, y0 in f[w, :2]:
            ph = np.arctan2(y0, x0); n_in = -np.array([np.cos(ph), np.sin(ph)]); tg = circ * np.array([-np.sin(ph), np.cos(ph)])   # along wall in circulation sense
            d0 = R - np.hypot(x0, y0); depths.append(d0)
            # local grid point (gx along, gy depth below the defect) -> lab position; the wall is at gy = -d0
            px = x0 + GX * tg[0] + GY * n_in[0]; py = y0 + GX * tg[1] + GY * n_in[1]
            ci = px + Q[0] - 0.5; cj = py + Q[1] - 0.5
            u = map_coordinates(U, [ci, cj], order=1, mode="nearest"); v = map_coordinates(V, [ci, cj], order=1, mode="nearest")
            inside = np.hypot(px, py) < R; ux = (u * tg[0] + v * tg[1]) * inside; uy = (u * n_in[0] + v * n_in[1]) * inside
            acc[0] += ux; acc[1] += uy; nacc += 1
    if nacc == 0: return None
    F = acc / nacc; d = np.median(depths)
    def grad_at(yd):
        j = int(np.argmin(np.abs(YS - yd))); i = int(np.argmin(np.abs(XS)))
        dux_dx = (F[0, i + 2, j] - F[0, i - 2, j]) / 2.0; duy_dx = (F[1, i + 2, j] - F[1, i - 2, j]) / 2.0
        dux_dy = (F[0, i, j + 2] - F[0, i, j - 2]) / 2.0; duy_dy = (F[1, i, j + 2] - F[1, i, j - 2]) / 2.0
        exx = 0.5 * (dux_dx - duy_dy); exy = 0.5 * (dux_dy + duy_dx); W = 0.5 * (dux_dy - duy_dx); E = np.hypot(exx, exy); phi = np.arctan2(exx, exy)
        tp = np.degrees(np.arccos(np.clip(W / (LAM * E), -1, 1)) - phi) if E > 0 else np.nan
        return dict(exx=exx, exy=exy, W=W, r=E / W if W else np.nan, phi=np.degrees(phi), tau_pred=tp, ux=F[0, i, j])
    prof = {yd: grad_at(yd) for yd in (1.0, 2.0, 3.0, 4.0, 5.0)}
    return dict(run="/".join(run.split("/")[-3:]), k=k, R=R, a=2 * np.pi * R / (2 * k - 2), depth=d,
                core=grad_at(0.0), prof=prof, tilt=np.degrees(np.median(tilts)) if len(tilts) > 5 else np.nan, F=F)


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/shapes_server/circle/circle/Req*/k*/s0*/DONE"))
    if os.environ.get("MAIN10"):   # 10-seed main circle sweep window (k 2-10, R 0.8-2.2 l_a) with fetched flow fields
        runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/flowfield.dat.gz"))
    TAG = os.environ.get("TAG", "")
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = [x for x in pool.map(one, runs) if x]
    res.sort(key=lambda x: x["a"])
    L = ["DEFECT-CENTRED WALL FLOW (shape-sweep circles, 3 seeds; local frame of each wall +1/2, x along the wall in the circulation sense, y into the cavity;",
         "averaged over wall defects and flow frames of the second half; gradients at the defect (centred, +-1 cell)). tau_pred = arccos(W/(lambda|E|)) - phi.",
         "run | k | spacing a | depth | at the defect: u_x, e_xy, e_xx, W, r = |E|/W, phi (deg), tau_pred | 2 cells deeper: r, phi, tau_pred | measured tilt"]
    for x in res:
        c, p2 = x["core"], x["prof"][2.0]
        L.append(f"{x['run']:14s} | {x['k']:2d} | {x['a']:5.1f} | {x['depth']:.2f} | {c['ux']:+.3f} {c['exy']:+.4f} {c['exx']:+.4f} {c['W']:+.4f} {c['r']:5.2f} {c['phi']:+6.1f} {c['tau_pred']:6.1f} | "
                 f"{p2['r']:5.2f} {p2['phi']:+6.1f} {p2['tau_pred']:6.1f} | {x['tilt']:5.1f}")
    a = np.array([x["a"] for x in res]); T = np.array([x["tilt"] for x in res]); TP = np.array([x["core"]["tau_pred"] for x in res])
    RR = np.array([x["core"]["r"] for x in res]); PH = np.array([x["core"]["phi"] for x in res]); ok = np.isfinite(T) & np.isfinite(TP)
    L.append(f"corr(tilt, log a) {np.corrcoef(T[ok], np.log(a[ok]))[0,1]:+.2f}; corr(tilt, r core) {np.corrcoef(T[ok], RR[ok])[0,1]:+.2f}; corr(tilt, phi core) {np.corrcoef(T[ok], PH[ok])[0,1]:+.2f}; "
             f"corr(tilt, tau_pred core) {np.corrcoef(T[ok], TP[ok])[0,1]:+.2f}; rms(tilt - tau_pred) {np.sqrt(np.mean((T[ok] - TP[ok]) ** 2)):.1f} deg")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.0))
    for x, col in ((res[0], "#c2452d"), (res[len(res) // 2], "0.4"), (res[-1], "#1f5fa8")):
        pass
    sel = [res[0], res[-1]]
    for j, (x, lab) in enumerate(zip(sel, (f"crowded, a = {res[0]['a']:.0f}", f"sparse, a = {res[-1]['a']:.0f}"))):
        F = x["F"]; sp = np.hypot(F[0], F[1]).max()
        ax[j].quiver(XS[::2], YS[::2], F[0][::2, ::2].T, F[1][::2, ::2].T, color="0.3", scale=8 * sp, width=0.005)
        ax[j].axhline(-x["depth"], color="k", lw=2); ax[j].plot(0, 0, "o", color="#c2452d"); ax[j].invert_yaxis(); ax[j].set_aspect("equal")
        ax[j].set_xlabel("along the wall, circulation sense (cells)"); ax[j].set_ylabel("depth below the defect (cells)"); ax[j].set_title(f"({'ab'[j]}) Flow around a wall +1/2\n{lab}", fontsize=9)
    sc = ax[2].scatter(a, RR, c=T, cmap="magma", vmin=40, vmax=70, s=18); ax[2].set_xscale("log"); ax[2].set_xlabel("wall-defect spacing (cells)"); ax[2].set_ylabel("strain / rotation at the defect, |E|/W")
    TP2 = np.array([x["prof"][2.0]["tau_pred"] for x in res]); PH2 = np.array([x["prof"][2.0]["phi"] for x in res])
    ax[2].cla(); sc = ax[2].scatter(a, PH2, c=T, cmap="magma", vmin=40, vmax=70, s=18); ax[2].set_xscale("log"); ax[2].set_xticks([7, 10, 15, 20, 30]); ax[2].set_xticklabels(["7", "10", "15", "20", "30"]); ax[2].minorticks_off()
    ax[2].set_xlabel("wall-defect spacing (cells)"); ax[2].set_ylabel("rotation of the strain axes, 2 cells inside (deg)"); ax[2].axhline(0, color="0.7", ls=":"); ax[2].set_title("(c) Local flow vs spacing", fontsize=10)
    ax[3].scatter(TP2, T, c=np.log(a), cmap="viridis", s=18); ax[3].plot([30, 90], [30, 90], "k--", lw=0.8); ax[3].set_xlabel("local Leslie tilt, 2 cells inside (deg)"); ax[3].set_ylabel("measured tilt (deg)")
    ax[3].set_title("(d) Local Leslie prediction", fontsize=10); fig.colorbar(sc, ax=ax[2], label="measured tilt (deg)")
    for a_ in ax[2:]: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_defect_flow{TAG}.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, f"defect_flow{TAG}.txt"), "w").write("\n".join(L) + "\n")
