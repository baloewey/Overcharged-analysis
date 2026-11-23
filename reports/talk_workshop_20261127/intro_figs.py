#!/usr/bin/env python3
"""Intro-slide graphics for the workshop talk (user, 23 Nov: softer intro for an active-matter audience that is not LC-expert).
  assets/rods.png       extensile active rods (force dipoles pushing outward along the axis)
  assets/defects.png    director field of +1/2 (comet, self-propelled towards its head in an extensile nematic) and -1/2 (trefoil, not self-propelled)
  assets/bulk_snap.png  MPCD active turbulence (runs/talk_intro/bulk80, act 0.07): director + defects, one frame
  assets/k1_snaps.png   k = 1 homeotropic circle (runs/talk_intro/k1_passive, k1_active): director + defects, passive vs active
Run: reports/figstyle/run.sh intro_figs.py [schematic|data]"""
import gzip, os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, "assets"); R = "/home/baloewe/Physics/Overcharged/runs/talk_intro"
RED, BLUE, GREY, ORANGE = "#c0392b", "#1f5fa8", "#6b6b6b", "#f08c00"


def sticks(ax, X, Y, TH, L, **kw):
    dx, dy = 0.5 * L * np.cos(TH), 0.5 * L * np.sin(TH)
    ax.plot(np.array([X - dx, X + dx]).reshape(2, -1), np.array([Y - dy, Y + dy]).reshape(2, -1), solid_capstyle="round", **kw)


def rods():
    fig, ax = plt.subplots(figsize=(4.2, 2.6)); rng = np.random.default_rng(3)
    for i in range(7):
        x0, y0, a = 1.7 * (i % 4) + 0.8 * (i // 4), 1.0 * (i // 4), 0.12 * rng.standard_normal()
        ux, uy = np.cos(a), np.sin(a); ax.plot([x0 - 0.45 * ux, x0 + 0.45 * ux], [y0 - 0.45 * uy, y0 + 0.45 * uy], "-", lw=9, color="#7fa7d6", solid_capstyle="round")
        if i == 5:
            for sgn in (1, -1):
                ax.annotate("", (x0 + sgn * 0.95 * ux, y0 + sgn * 0.95 * uy), (x0 + sgn * 0.5 * ux, y0 + sgn * 0.5 * uy),
                            arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.5, mutation_scale=16))
    ax.set_xlim(-0.6, 6.4); ax.set_ylim(-0.5, 1.5); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(os.path.join(A, "rods.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)


def defects():
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.9)); g = np.linspace(-1, 1, 15); X, Y = np.meshgrid(g, g)
    for ax, q, th0, title, col in ((axs[0], 0.5, 0.0, r"$+\frac{1}{2}$: self-propelled", RED), (axs[1], -0.5, 0.75 * np.pi, r"$-\frac{1}{2}$: not self-propelled", BLUE)):
        r_ = np.hypot(X, Y); m = (r_ > 0.12) & (r_ < 1.08); TH = q * np.arctan2(Y, X) + th0      # th0: -1/2 spokes at 90, 210, 330 deg
        sticks(ax, X[m], Y[m], TH[m], 0.13, color=GREY, lw=2.2)
        ax.plot(0, 0, "o", ms=13, color=col, mec="white", mew=1.5, zorder=5)
        if q > 0: ax.annotate("", (-0.72, 0), (-0.12, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=4, mutation_scale=26), zorder=6)
        ax.set_title(title, fontsize=15, color=col); ax.set_xlim(-1.12, 1.12); ax.set_ylim(-1.12, 1.12); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(os.path.join(A, "defects.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)


def read_all(run):
    """all director frames {t: (x, y, nx, ny, S)} and defect frames {t: (x, y, q, angle)}"""
    def opn(name):
        fn = os.path.join(R, run, name + ".gz"); fn = fn if os.path.exists(fn) else fn[:-3]; return (gzip.open if fn.endswith(".gz") else open)(fn, "rt")
    D, cur, tcur = {}, [], None
    with opn("directorfield.dat") as f:
        for line in f:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            t = float(p[0])
            if t != tcur:
                if cur: D[tcur] = np.array(cur)
                cur, tcur = [], t
            cur.append((float(p[1]), float(p[2]), float(p[4]), float(p[5]), float(p[7])))
    if cur: D[tcur] = np.array(cur)
    F, cur, tcur, nexp = {}, [], None, 0
    with opn("defects.dat") as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                try: t, n = float(p[0]), int(p[1])
                except ValueError: continue
                if tcur is not None and len(cur) == nexp: F[tcur] = np.array(cur).reshape(-1, 4)
                cur, tcur, nexp = [], t, n
            elif len(p) == 4 and tcur is not None: cur.append([float(v) for v in p])
    if tcur is not None and len(cur) == nexp: F[tcur] = np.array(cur).reshape(-1, 4)
    return D, F


def draw_field(ax, Dt, Ft, circle=None, step=2, L=0.8):
    x, y, nx, ny, S = Dt.T; m = S > 0.05
    if circle is not None:
        cx, cy, Rr = circle; m &= np.hypot(x + 0.5 - cx, y + 0.5 - cy) < Rr - 0.5
    m &= (x % step == 0) & (y % step == 0)
    sticks(ax, x[m] + 0.5, y[m] + 0.5, np.arctan2(ny[m], nx[m]), L * step, color="#9aa3ad", lw=1.3)
    for xd, yd, q, a in Ft:
        if q > 0:
            ax.plot([xd, xd + 2.2 * np.cos(a)], [yd, yd + 2.2 * np.sin(a)], "-", color=RED, lw=4, solid_capstyle="round", zorder=4)
            ax.plot(xd, yd, "o", ms=9, color=RED, mec="white", mew=1.3, zorder=5)
        else:
            for a_ in (0, 2 * np.pi / 3, -2 * np.pi / 3):
                ax.plot([xd, xd + 1.6 * np.cos(a + a_)], [yd, yd + 1.6 * np.sin(a + a_)], "-", color=BLUE, lw=3.2, solid_capstyle="round", zorder=4)
    if circle is not None: ax.add_patch(plt.Circle(circle[:2], circle[2], fill=False, lw=3, color="k"))
    ax.set_aspect("equal"); ax.axis("off")


def frame(D, F, t, box, circle=None, size=5.0, title=None):
    fig = plt.figure(figsize=(size, size * 1.08 if title else size)); ax = fig.add_axes([0, 0, 1, 0.92 if title else 1])
    draw_field(ax, D[t], F.get(t, np.zeros((0, 4))), circle=circle); ax.set_xlim(*box[:2]); ax.set_ylim(*box[2:])
    if title: fig.text(0.5, 0.95, title, ha="center", va="center", fontsize=17)
    return fig


def movie(D, F, name, box, circle=None, nmax=150, fps=10):
    import subprocess
    fd = os.path.join(A, name + "_frames"); os.makedirs(fd, exist_ok=True); [os.remove(os.path.join(fd, f)) for f in os.listdir(fd)]
    ts = sorted(t for t in D if t in F)[-nmax:]
    for i, t in enumerate(ts):
        fig = frame(D, F, t, box, circle); fig.savefig(os.path.join(fd, f"f{i:03d}.png"), dpi=110, facecolor="white"); plt.close(fig)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", os.path.join(fd, "f%03d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    "-pix_fmt", "yuv420p", "-c:v", "libx264", os.path.join(A, name + ".mp4")], check=True)
    subprocess.run(["cp", os.path.join(fd, "f000.png"), os.path.join(A, name + "_poster.png")], check=True); print("movie", name, len(ts), "frames", ts[0], "-", ts[-1])


def data():
    import json
    D, F = read_all("bulk80"); movie(D, F, "bulk_movie", (0, 80, 0, 80))
    bc = json.load(open(os.path.join(R, "k1_active", "input.json")))["BC"][0]; circ = (bc["Q"][0], bc["Q"][1], bc["R"]); box = (circ[0] - circ[2] - 2, circ[0] + circ[2] + 2, circ[1] - circ[2] - 2, circ[1] + circ[2] + 2)
    Dp, Fp = read_all("k1_passive"); tp = max(t for t in Dp if t in Fp); fig = frame(Dp, Fp, tp, box, circ, size=4.4); fig.savefig(os.path.join(A, "k1_passive.png"), dpi=250, facecolor="white"); plt.close(fig)
    Da, Fa = read_all("k1_active"); movie(Da, Fa, "k1_active_movie", box, circ)
    S = {}
    for lab, Ff in (("passive", Fp), ("active", Fa)):           # defect census: wall (< 4 cells from the wall) vs bulk, +1/2 and -1/2
        rows = []
        for t, f in Ff.items():
            if not len(f): rows.append((0, 0, 0)); continue
            d = circ[2] - np.hypot(f[:, 0] - circ[0], f[:, 1] - circ[1]); rows.append((((d < 4) & (f[:, 2] > 0)).sum(), ((d >= 4) & (f[:, 2] > 0)).sum(), (f[:, 2] < 0).sum()))
        rows = np.array(rows[len(rows) // 5:]); S[lab] = rows.mean(0)
        print(f"k1 {lab}: frames {len(rows)} | wall +1/2 {rows[:, 0].mean():.2f} | bulk +1/2 {rows[:, 1].mean():.2f} | -1/2 {rows[:, 2].mean():.2f} | net bulk charge {(0.5 * rows[:, 1] - 0.5 * rows[:, 2]).mean():.2f}")
    open(os.path.join(A, "k1_census.txt"), "w").write("\n".join(f"{k}: wall+ {v[0]:.2f} bulk+ {v[1]:.2f} minus {v[2]:.2f}" for k, v in S.items()) + "\n")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "schematic"
    if what == "schematic": rods(); defects()
    else: data()
