#!/usr/bin/env python3
"""Movie: crowded wall layers tilt less in MPCD, not in the continuum (user, 2026-10-16). 2 x 2 cavities, R 31.83, lambda 2:
top = active MPCD (k 4, spacing 33 cells; k 16, spacing 6.7 cells), 300 t.u. of defect frames (every 2 t.u.) in a coherent stretch of the second half;
bottom = active continuum Q-tensor + backflow (one constant, zeta 0.3; same k), 300 t.u. after 3000 (shear_tilt/active_snapshots.py), director sticks and |S|.
Each wall +1/2 is drawn as a comet (dot = core, tail along its radial arm) and a grey tick along the local outward normal; colour = its tilt (0-90 deg).
Bottom strip: layer tilt |arg<exp(i tau)>| against frame for the four cases, dashed arccos(1/2) = 60 deg. -> movie_tilt_crowding.mp4"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); ROOT = "/home/baloewe/Physics/Overcharged"
sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/shear_tilt")
from edge_criteria import read_with_angle  # noqa
import shear_tilt as M  # noqa
NF = 150; cmap = plt.cm.plasma


def mpcd(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
    C = []
    for f in fr:
        r = np.hypot(f[:, 0], f[:, 1]) if len(f) else np.zeros(0); w = (f[:, 2] > 0) & (R - r < 4) if len(f) else np.zeros(0, bool)
        C.append(abs(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))) if w.sum() >= 3 else 0)
    C = np.array(C); s = int(np.argmax(np.convolve(C > 0.9, np.ones(NF), "valid")))     # most coherent 300 t.u. stretch
    frames = []
    for f in fr[s:s + NF]:
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4); frames.append((f, w))
    return dict(kind="mpcd", R=R, k=k, frames=frames)


def cont(k):
    Z = np.load(os.path.join(ROOT, "reports/master/theory/shear_tilt", f"snaps_k{k:02d}.npz")); X, Y, R = Z["X"], Z["Y"], float(Z["R"]); frames = []
    for n in range(min(NF, len(Z["t"]))):
        q1, q2 = Z["q1"][n].astype(float), Z["q2"][n].astype(float); D = M.defects(q1, q2)
        f = np.array([(px, py, q, a) for px, py, q, a in D]) if D else np.zeros((0, 4)); r = np.hypot(f[:, 0], f[:, 1]) if len(f) else np.zeros(0)
        w = (f[:, 2] > 0) & (r > R - 6) if len(f) else np.zeros(0, bool); frames.append((f, w, q1, q2))
    return dict(kind="cont", R=R, k=k, frames=frames, X=X, Y=Y)


def layer_tilt(f, w):
    if w.sum() < 3: return np.nan
    return np.degrees(abs(np.angle(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))))))


if __name__ == "__main__":
    S = [mpcd(ROOT + "/runs/cavity/server_main/main_act0.07/rho0.04/k04/s01"), mpcd(ROOT + "/runs/cavity/server_main/main_act0.07/rho0.16/k16/s01"), cont(4), cont(16)]
    titles = ["MPCD, k = 4 (6 wall defects, spacing 33 cells)", "MPCD, k = 16 (30 wall defects, spacing 6.7 cells)",
              "Continuum, k = 4 (one constant, active, backflow)", "Continuum, k = 16 (one constant, active, backflow)"]
    TS = [[layer_tilt(*fr[:2]) for fr in s["frames"]] for s in S]
    out = os.path.join(H, "frames_tiltmovie"); os.makedirs(out, exist_ok=True)
    for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
    for n in range(NF):
        fig = plt.figure(figsize=(11, 12.5)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.55], hspace=0.28, wspace=0.08)
        for i, (s, tt) in enumerate(zip(S, titles)):
            ax = fig.add_subplot(gs[i // 2, i % 2]); R = s["R"]; ph = np.linspace(0, 2 * np.pi, 300); ax.plot(R * np.cos(ph), R * np.sin(ph), "k-", lw=1.5)
            fr = s["frames"][min(n, len(s["frames"]) - 1)]; f, w = fr[0], fr[1]
            if s["kind"] == "cont":
                q1, q2 = fr[2], fr[3]; X, Y = s["X"], s["Y"]; inside = np.hypot(X, Y) < R; Sm = np.where(inside, 2 * np.hypot(q1, q2), np.nan); th = 0.5 * np.arctan2(q2, q1)
                ax.pcolormesh(X, Y, Sm, cmap="Greys_r", vmin=0, vmax=0.9, shading="nearest", rasterized=True)
                st = 5; m = inside[::st, ::st]; xs, ys, tq = X[::st, ::st][m], Y[::st, ::st][m], th[::st, ::st][m]
                ax.plot(np.c_[xs - 1.3 * np.cos(tq), xs + 1.3 * np.cos(tq)].T, np.c_[ys - 1.3 * np.sin(tq), ys + 1.3 * np.sin(tq)].T, color="#8fb3d9", lw=0.6)
            for (x, y, q, a), iw in zip(f, w):
                if q > 0:
                    phi = np.arctan2(y, x); tau = np.degrees(abs(np.angle(np.exp(1j * (a - phi))))) if iw else None
                    col = cmap(min(tau, 90) / 90) if iw else "0.5"
                    if iw: ax.plot([x, x + 4 * np.cos(phi)], [y, y + 4 * np.sin(phi)], color="0.55", lw=0.8)
                    ax.plot([x, x + 4 * np.cos(a)], [y, y + 4 * np.sin(a)], color=col, lw=2.2); ax.plot(x, y, "o", color=col, ms=5, mec="k", mew=0.4)
                elif q < 0:
                    for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 2 * np.cos(b)], [y, y + 2 * np.sin(b)], color="#1f5fa8", lw=1.2)
            tl = TS[i][min(n, len(TS[i]) - 1)]; ax.set_title(tt + (f"\nlayer tilt {tl:4.1f}°" if np.isfinite(tl) else "\nlayer tilt –"), fontsize=10)
            ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-R - 5, R + 5); ax.set_ylim(-R - 5, R + 5)
        a2 = fig.add_subplot(gs[2, :])
        for i, (c, ls, lab) in enumerate(((plt.cm.plasma(0.15), "-", "MPCD k 4"), (plt.cm.plasma(0.15), "--", "MPCD k 16"), ("0.2", "-", "continuum k 4"), ("0.2", "--", "continuum k 16"))):
            a2.plot(np.arange(len(TS[i])) * 2, TS[i], ls, color=c, lw=1.5, label=lab)
        a2.axhline(60, color="k", ls=":", lw=0.8); a2.axvline(2 * n, color="#c2452d", lw=1); a2.set_ylim(30, 80); a2.set_xlim(0, 2 * NF)
        a2.set_xlabel("time in the window (t.u.; MPCD and continuum units)"); a2.set_ylabel("layer tilt (deg)"); a2.legend(fontsize=8, frameon=False, ncol=4, loc="lower center")
        [a2.spines[s_].set_visible(False) for s_ in ("top", "right")]
        cax = fig.add_axes([0.92, 0.42, 0.012, 0.25]); fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 90)), cax=cax, label="wall-defect tilt (deg)")
        fig.suptitle("Crowded wall layers: MPCD tilt falls, continuum tilt stays at arccos(1/λ)", fontsize=12, y=0.93)
        fig.savefig(os.path.join(out, f"f{n:04d}.png"), dpi=90); plt.close(fig)
    mp4 = os.path.join(H, "movie_tilt_crowding.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", os.path.join(out, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, [np.nanmedian(t) for t in TS])
