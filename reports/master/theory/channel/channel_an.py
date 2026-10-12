#!/usr/bin/env python3
"""Analysis of winding-anchored periodic channels (branch channelwind, eMode 4). Per run: defects per frame (cavity coordinates = lattice coordinates), wall distance
d = min(y, W - y); +1/2 within 4 cells of the bottom / top wall, bulk charge (|d| >= 4), total charge; tilt of wall +1/2 relative to the wall normal; from the flow frames:
vorticity along the channel mid-line, number of sign changes (vortex count) and dominant wavenumber; mean shear profile <v_x>(y). Optional snapshot figure.
Usage: channel_an.py RUNDIR [...] [--fig]"""
import os, sys, json, gzip, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
def opn(fn): return gzip.open(fn + ".gz", "rt") if os.path.exists(fn + ".gz") else open(fn)
def flow_frames(run):
    rows = np.loadtxt((l for l in opn(run + "/flowfield.dat") if len(l.split()) == 7 and l.split()[0][0].isdigit())); return rows
for run in [a for a in sys.argv[1:] if not a.startswith("--")]:
    J = json.load(open(run + "/input.json")); L, W = J["domain"][:2]; m = J["BC"][2].get("customAnch", 0)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); h = len(fr) // 2
    nb, nt, qb, qt, tilts = [], [], [], [], []
    for f in fr[h:]:
        if not len(f): continue
        y = f[:, 1]; d = np.minimum(y, W - y); bot = (f[:, 2] > 0) & (y < 4); top = (f[:, 2] > 0) & (W - y < 4); bulk = d >= 4
        nb.append(bot.sum()); nt.append(top.sum()); qb.append(0.5 * np.sign(f[bulk, 2]).sum()); qt.append(0.5 * np.sign(f[:, 2]).sum())
        for i in np.flatnonzero(bot): tilts.append(np.degrees(abs(np.angle(np.exp(1j * (f[i, 3] - np.pi / 2))))))
        for i in np.flatnonzero(top): tilts.append(np.degrees(abs(np.angle(np.exp(1j * (f[i, 3] + np.pi / 2))))))
    line = f"{os.path.basename(run)}: L {L} W {W} m {m} (charge 2m = {2 * m}); t {t[h]:.0f}-{t[-1]:.0f}: wall +1/2 bottom {np.mean(nb):.2f} top {np.mean(nt):.2f} (2m per wall if all at the walls: {2 * m}); bulk charge {np.mean(qb):+.2f}; total {np.mean(qt):+.2f}; wall-defect |tilt| median {np.median(tilts) if tilts else np.nan:.1f} deg"
    try:
        R = flow_frames(run); ft = np.unique(R[:, 0]); ft = ft[ft >= t[h]]; ncross, kdom, prof = [], [], []
        for tt in ft:
            r = R[R[:, 0] == tt]; ix, iy = r[:, 1].astype(int), r[:, 2].astype(int); U = np.zeros((L, W)); V = np.zeros((L, W)); U[ix, iy] = r[:, 4]; V[ix, iy] = r[:, 5]
            om = np.gradient(V, axis=0) - np.gradient(U, axis=1); mid = om[:, W // 4: 3 * W // 4].mean(1); mid -= mid.mean()
            ncross.append(np.sum(np.sign(mid) != np.sign(np.roll(mid, 1)))); spec = np.abs(np.fft.rfft(mid)) ** 2; kdom.append(np.argmax(spec[1:]) + 1); prof.append(U.mean(0))
        line += f"; mid-line vorticity sign changes {np.mean(ncross):.1f} (vortices ~ {np.mean(ncross):.1f}), dominant wavenumber {np.median(kdom):.0f} per L; |<v_x>| near walls {np.abs(np.mean(prof, 0)[[1, W - 2]]).mean():.3f}, centre {np.mean(prof, 0)[W // 2]:+.3f}"
    except Exception as ex: line += f"; flow: {ex}"
    print(line)
    if "--fig" in sys.argv:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        f = fr[-1]; fig, ax = plt.subplots(figsize=(10, 10 * W / L + 0.6))
        try:
            D = np.loadtxt((l for l in opn(run + "/directorfield.dat") if len(l.split()) >= 7 and l.split()[0][0].isdigit())); tl = D[:, 0].max(); D = D[D[:, 0] == tl]
            ax.quiver(D[:, 1] + 0.5, D[:, 2] + 0.5, D[:, 4], D[:, 5], headwidth=0, headlength=0, headaxislength=0, pivot="mid", scale=70, width=0.0015, color="0.4")
        except Exception as ex: print("director:", ex)
        for x, y, q, a in f:
            if q > 0: ax.plot(x, y, "o", color="#c2452d", ms=6); ax.annotate("", (x - 2 * np.cos(a), y - 2 * np.sin(a)), (x, y), arrowprops=dict(arrowstyle="-|>", color="#c2452d"))
            else: ax.plot(x, y, "^", color="#1f5fa8", ms=6)
        ax.axhline(0, color="k", lw=3); ax.axhline(W, color="k", lw=3); ax.set_xlim(0, L); ax.set_ylim(-1, W + 1); ax.set_aspect("equal"); ax.set_title(os.path.basename(run) + f"  t = {t[-1]:.0f}")
        fig.savefig(run + "/snapshot.png", dpi=110, bbox_inches="tight"); plt.close(fig)
