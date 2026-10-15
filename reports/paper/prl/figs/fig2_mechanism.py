#!/usr/bin/env python3
"""PRL Fig. 2 (mechanism): flow-induced ejection. Single column. Run: figs/run.sh fig2_mechanism.py
(a) Elastic energy (k = 6) above the 2k ring with its natural, winding tilts, when every wall defect is forced to a common tilt tau: the 2k ring itself
    ("ring frustration") and the ejected structure (2k-2 ring + bulk pair; = frustration + E_B - E_A). Source master/theory/ejection/ejection_energetics.txt;
    line = anisotropic K1/K3 = 1.82, core 1.5; band = both elasticities and cores 1.5/2.5. At tau = 0 the ejected value is the ejection cost without flow.
(c) Passive MPCD cavity (k = 6) under an imposed vortex of amplitude A: wall +1/2 count and tilt coherence C vs time (runs/swirl_threshold, swirl_tumble;
    same reduction as theory_report/make_figs.py: wall = within 4 cells, every 5th frame, 10-point running mean).
(b) Active continuum Q-tensor model with backflow (no imposed flow; zeta 0.3, lambda 2, k 6), final state from the passive 12-ring
    (master/theory/shear_tilt/state2_active_lam2.0_z0.3_A0.0.npz): order S, director, flow, defects (shear_tilt.defects)."""
import os, sys, json, re
import numpy as np
from scipy.spatial import cKDTree
from prlstyle import plt, COL, RED, BLUE, GREEN, GREY, DARK, ROOT, label, save
sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/shear_tilt")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
import shear_tilt as ST  # noqa

# ---------- (a) data ----------
def parse_ts():
    out = {}
    for line in open(ROOT + "/reports/master/theory/ejection/tilted_structures.txt"):
        if "Delta E = E_B - E_A =" not in line: continue
        kind = "aniso" if line.startswith("anisotropic") else "one"; s = float(line.split("s =")[1].split("R")[0])
        seg = line.split("Delta E = E_B - E_A =")[1].split("kBT")[0].split()
        out[(kind, s)] = np.array([[float(v) for v in p.split(":")] for p in seg])
    return out

TAUS = np.array([0, 30, 60, 90])
def parse_ej(k):
    """ejection_energetics.txt: per (elasticity, core) the frustration of the 2k ring forced to a common tilt tau (relative to its natural tilts) and
    E_B - E_A; the ejected structure relative to the natural ring = frustration + (E_B - E_A)."""
    out = {}
    for line in open(ROOT + "/reports/master/theory/ejection/ejection_energetics.txt"):
        m = re.match(r"k=(\d+) (\S+)\s+core ([\d.]+):", line)
        if not m or int(m.group(1)) != k: continue
        fr = [float(x) for x in re.findall(r"ring frustration ([+-]\d+)", line)]; ba = [float(x) for x in re.findall(r"B-A ([+-]\d+)", line)]
        out[(m.group(2), float(m.group(3)))] = dict(ring=np.array(fr), eject=np.array(fr) + np.array(ba))
    return out

# ---------- (b) data ----------
def series(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); T, NW, CO = [], [], []
    for i in range(0, len(t), 5):
        f = fr[i]
        if not len(f): continue
        d, j = tree.query(f[:, :2]); w = (d < 4) & (f[:, 2] > 0)
        T.append(t[i]); NW.append(w.sum()); CO.append(abs(np.mean(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]]))))) if w.any() else np.nan)
    sm = lambda a: np.convolve(np.nan_to_num(a), np.ones(10) / 10, mode="same")
    return np.array(T)[5:-5], sm(np.array(NW, float))[5:-5], sm(np.array(CO, float))[5:-5]

RUNS = [("swirl_threshold/passive_A0.001_s1", "$A=0.001$", "0.65", "-"), ("swirl_threshold/passive_A0.003_s2", "$A=0.003$", "0.35", "-"),
        ("swirl_threshold/passive_A0.01_s2", "$A=0.01$", BLUE, "-"), ("swirl_tumble/passive_lam4.0_A0.01", "$A=0.01$, $\\lambda=4$", GREEN, "-"),
        ("swirl_tumble/passive_lam0.5_A0.01", "$A=0.01$, $\\lambda=0.5$ (tumbling)", RED, "-")]

if __name__ == "__main__":
    fig = plt.figure(figsize=(COL, 4.75))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.35, 0.36, 0.8, 0.8], width_ratios=[1.0, 0.95], hspace=0.12, wspace=0.42)
    # (a) energies relative to the 2k ring with its natural (winding) tilts, k = 6
    ax = fig.add_subplot(gs[0, 0]); D = parse_ej(6)
    for key, col, lab in (("ring", DARK, "$2k$ ring"), ("eject", RED, "Ejected ($2k-2$ + pair)")):
        A = np.array([D[v][key] for v in D]); ref = D[("anisotropic", 1.5)][key]
        ax.fill_between(TAUS, A.min(0), A.max(0), color=col, alpha=0.18, lw=0)
        ax.plot(TAUS, ref, "-o", color=col, ms=2.5, label=lab)
    ax.set_xlim(0, 90); ax.set_xticks([0, 30, 60, 90]); ax.set_ylim(0, 540); ax.set_yticks([0, 100, 200, 300, 400, 500])
    ax.set_xlabel(r"Imposed common tilt $\tau$ (deg)", labelpad=1); ax.set_ylabel("Energy above the natural\n$2k$ ring ($k_BT$)", labelpad=2)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=6, handlelength=1.4); label(ax, "a", x=-0.3)
    # (c)
    ac = fig.add_subplot(gs[0, 1]); Z = np.load(ROOT + "/reports/master/theory/shear_tilt/state2_active_lam2.0_z0.3_A0.0.npz")
    q1, q2, vx, vy = Z["q1"], Z["q2"], Z["vx"], Z["vy"]; X, Y, R, inside = ST.X, ST.Y, ST.R, ST.inside
    S = np.where(inside, 2 * np.hypot(q1, q2), np.nan); th = 0.5 * np.arctan2(q2, q1)
    ac.pcolormesh(X, Y, S, cmap="Greys_r", vmin=0, vmax=ST.S0, shading="nearest", rasterized=True)
    st = 6; m = inside[::st, ::st] & (np.hypot(X, Y)[::st, ::st] < R - 1.5); L = 1.6
    xs, ys, tt = X[::st, ::st][m], Y[::st, ::st][m], th[::st, ::st][m]
    ac.plot(np.c_[xs - L * np.cos(tt), xs + L * np.cos(tt)].T, np.c_[ys - L * np.sin(tt), ys + L * np.sin(tt)].T, color="0.55", lw=0.45)
    sv = 8; mv = inside[::sv, ::sv] & (np.hypot(X, Y)[::sv, ::sv] < R - 3)
    ac.quiver(X[::sv, ::sv][mv], Y[::sv, ::sv][mv], vx[::sv, ::sv][mv], vy[::sv, ::sv][mv], color="#d9962b", scale=9, width=0.008, headwidth=3.5, headlength=3.5, zorder=3)
    ph = np.linspace(0, 2 * np.pi, 300); ac.plot(R * np.cos(ph), R * np.sin(ph), color="k", lw=1)
    nw = nb = 0
    for px, py, q, ang in ST.defects(q1, q2):
        if q > 0:
            ac.plot(px, py, "o", color=RED, ms=3.2, mec="w", mew=0.4, zorder=5); ac.plot([px, px - 4 * np.cos(ang)], [py, py - 4 * np.sin(ang)], color=RED, lw=1.1, zorder=4)
            nw += np.hypot(px, py) > R - 6; nb += np.hypot(px, py) <= R - 6
        elif q < 0: ac.plot(px, py, "^", color=BLUE, ms=3.2, zorder=5)
    ac.set_aspect("equal"); ac.axis("off"); ac.set_xlim(-R - 1, R + 1); ac.set_ylim(-R - 1, R + 1)
    ac.text(0, -R - 3, f"Active continuum model:\n{nw} wall + {nb} bulk $+\\frac{{1}}{{2}}$", ha="center", va="top", fontsize=6.5); label(ac, "b", x=0.02, y=0.97)
    print("continuum snapshot: wall", nw, "bulk", nb)
    # (b)
    a1 = fig.add_subplot(gs[2, :]); a2 = fig.add_subplot(gs[3, :], sharex=a1)
    for rel, lab, col, ls in RUNS:
        run = os.path.join(ROOT, "runs", rel); T, NW, CO = series(run)
        a1.plot(T, NW, ls, color=col, lw=0.9, label=lab); a2.plot(T, CO, ls, color=col, lw=0.9)
    for v in (10, 12): a1.axhline(v, color="0.75", ls=":", lw=0.6)
    a1.set_yticks([10, 11, 12]); a1.set_ylim(9.6, 12.4); a1.set_ylabel("Wall $+\\frac{1}{2}$", labelpad=2)
    a1.text(3010, 12, "$2k$", fontsize=6.5, va="center"); a1.text(3010, 10, "$2k-2$", fontsize=6.5, va="center")
    plt.setp(a1.get_xticklabels(), visible=False); a2.set_ylim(0, 1.05); a2.set_ylabel("Tilt coherence $C$", labelpad=2); a2.set_xlabel("Time (t.u.)", labelpad=1)
    a1.set_xlim(0, 3000)
    h, l = a1.get_legend_handles_labels()
    a2.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=3, fontsize=6, columnspacing=1.0, handlelength=1.4,
              title="Passive MPCD, $k=6$, imposed vortex of amplitude $A$ ($\\lambda=2$ unless noted)", title_fontsize=6)
    label(a1, "c", x=-0.1, y=0.98)
    save(fig, "fig2_mechanism")
