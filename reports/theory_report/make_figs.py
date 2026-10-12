#!/usr/bin/env python3
"""Figures for the theory report (new ones only; the rest are copied from reports/master/figs).
fig_lock.png: symmetry/lock schematic; fig_ejection.png: passive cavities under imposed swirl; fig_shapes_mech.png: escape, coherence, cascade tests."""
import csv, glob, json, os, pickle, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged"; sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/avalanche")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); F = os.path.join(H, "figs")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
C1, C2, C3 = "#1f5fa8", "#c2452d", "#3a8a4f"
# ---------------- 1. lock schematic, k = 4 ----------------
def lock():
    k = 4; fig, ax = plt.subplots(1, 3, figsize=(10, 4.0))
    for a in ax: a.set_aspect("equal"); a.axis("off"); a.set_xlim(-1.25, 1.25); a.set_ylim(-1.25, 1.25)
    ph = np.linspace(0, 2 * np.pi, 400); 
    for a in ax: a.plot(np.cos(ph), np.sin(ph), color="0.2", lw=1.2)
    for p in np.linspace(0, 2 * np.pi, 48, endpoint=False):           # wall director theta = k phi
        th = k * p; c = np.array([np.cos(p), np.sin(p)]) * 1.08; d = 0.07 * np.array([np.cos(th), np.sin(th)])
        for a in ax: a.plot([c[0] - d[0], c[0] + d[0]], [c[1] - d[1], c[1] + d[1]], color="0.45", lw=1.3)
    hp = 2 * np.pi * np.arange(2 * k - 2) / (2 * k - 2)
    ax[0].scatter(1.08 * np.cos(hp), 1.08 * np.sin(hp), s=55, facecolor="none", edgecolor=C2, lw=1.5, zorder=5)
    ax[0].set_title(f"(a) anchoring $\\theta=k\\phi$, $k={k}$:\n{2*k-2} homeotropic points (circles)", fontsize=9)
    def defects(a, phis, tilts, col):
        r = 0.8
        for p, t in zip(phis, tilts):
            x, y = r * np.cos(p), r * np.sin(p); ang = p + t
            a.plot(x, y, "o", color=col, ms=5)
            a.annotate("", xy=(x + 0.2 * np.cos(ang), y + 0.2 * np.sin(ang)), xytext=(x, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3))
    defects(ax[1], hp, np.zeros_like(hp), C1)
    ax[1].set_title(f"(b) active: $2k-2={2*k-2}$ $+\\frac{{1}}{{2}}$ on the sites,\ntails (arrows) at the wall, tilt 0 (commensurate)", fontsize=9)
    p8 = 2 * np.pi * (np.arange(2 * k) + 0.5) / (2 * k); t8 = np.angle(np.exp(1j * (2 * (k - 1) * p8 + 0.4)))
    defects(ax[2], p8, t8, C2)
    ax[2].set_title(f"(c) passive: $2k={2*k}$ $+\\frac{{1}}{{2}}$; tilt $=2(k-1)\\phi_j$:\nsteps of $-2\\pi/k$, cannot all vanish", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(os.path.join(F, "fig_lock.png"), dpi=200); plt.close(fig)
# ---------------- 2. passive ejection under imposed swirl ----------------
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
    return np.array(T), sm(np.array(NW, float)), sm(np.array(CO, float))
def ejection():
    runs = [("swirl_threshold/passive_A0.001_s1", "swirl $A=0.001$, $\\lambda=2$", "0.6"), ("swirl_threshold/passive_A0.003_s2", "$A=0.003$, $\\lambda=2$", "0.35"),
            ("swirl_threshold/passive_A0.01_s2", "$A=0.01$, $\\lambda=2$", C1), ("swirl_tumble/passive_lam4.0_A0.01", "$A=0.01$, $\\lambda=4$", C3),
            ("swirl_tumble/passive_lam0.5_A0.01", "$A=0.01$, $\\lambda=0.5$ (tumbling)", C2)]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.3))
    for rel, lab, col in runs:
        run = os.path.join(ROOT, "runs", rel)
        if not os.path.exists(run + "/input.json"): continue
        T, NW, CO = series(run); m = slice(5, -5)
        ax[0].plot(T[m], NW[m], color=col, lw=1.3, label=lab); ax[1].plot(T[m], CO[m], color=col, lw=1.3)
    ax[0].axhline(12, color="0.7", ls=":", lw=0.8); ax[0].axhline(10, color="0.7", ls=":", lw=0.8)
    ax[0].set_ylabel("wall $+\\frac{1}{2}$ (passive, $k=6$)"); ax[1].set_ylabel("tilt coherence of the wall $+\\frac{1}{2}$")
    for a in ax: a.set_xlabel("time")
    fig.legend(*ax[0].get_legend_handles_labels(), fontsize=7.5, frameon=False, loc="lower center", ncol=5); ax[0].set_title("(a) wall count: $2k=12\\to2k-2=10$", fontsize=9)
    ax[1].set_title("(b) ejection follows a common tilt", fontsize=9)
    fig.tight_layout(rect=[0, 0.08, 1, 1]); fig.savefig(os.path.join(F, "fig_ejection.png"), dpi=200); plt.close(fig)
# ---------------- 3. shape mechanism ----------------
def shapes():
    S = ROOT + "/reports/master/theory"
    rows = list(csv.DictReader(open(S + "/shapes/kinetics_carriers.csv")))
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4)); cols = {1.0: C1, 1.5: "#d08a1e", 2.0: C3}
    for R in (1.0, 1.5, 2.0):
        for fam, ls in (("ellipse", "-"), ("oval", "--")):
            sel = [r for r in rows if r["family"] in (fam, "circle") and float(r["Req_la"]) == R]
            asp = sorted(set(float(r["aspect"]) for r in sel))
            g = [np.mean([float(r["r_esc"]) / float(r["nwall"]) for r in sel if float(r["aspect"]) == a]) for a in asp]
            ax[0].plot(asp, g, ls, marker="o", ms=3, color=cols[R], label=f"$R_\\mathrm{{eq}}={R}\\,\\ell_a$" if fam == "ellipse" else None)
    ax[0].set_xlabel("aspect ratio"); ax[0].set_ylabel("escapes per wall $+\\frac{1}{2}$ per 1000 t.u."); ax[0].set_title("(a) escape rate (solid ellipses, dashed ovals)", fontsize=9)
    ax[0].legend(fontsize=7, frameon=False)
    D = pickle.load(open(S + "/shapes/interval_coupling_cond.pkl", "rb")); man = {r["run"]: r for r in csv.DictReader(open(S + "/shapes/shape_runs.csv"))}
    bins = [(0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01)]; xc = [0.25, 0.6, 0.775, 0.93]
    for key, mk, lab in (("RE", "o-", "elongated"), ("RC", "s--", "circle-like")):
        for R in (1.0, 1.5, 2.0):
            RS = [r for r in D[key] if "error" not in r and len(r["NW"]) > 10 and float(man[r["rel"]]["Req_la"]) == R]
            co = np.concatenate([r["CO"] for r in RS]); vac = np.concatenate([(2 * r["k"] - 2) - r["NW"] for r in RS])
            y = [np.nanmean(vac[(co >= a) & (co < b)]) if ((co >= a) & (co < b)).sum() > 30 else np.nan for a, b in bins]
            ax[1].plot(xc, y, mk, ms=3.5, color=cols[R], label=f"{lab}" if R == 2.0 else None, mfc="white" if key == "RC" else cols[R])
    ax[1].set_xlabel("tilt coherence in the 100 t.u. interval"); ax[1].set_ylabel("vacancies"); ax[1].legend(fontsize=7, frameon=False)
    ax[1].set_title("(b) transient coherence does not refill", fontsize=9)
    import avalanche as AV
    res = [r for r in pickle.load(open(S + "/avalanche/avalanche_events.pkl", "rb")) if "error" not in r]
    mm = {r["dir"]: r for r in csv.DictReader(open(ROOT + "/runs/shapes_server/manifest_shape_sweep.tsv"), delimiter="\t")}
    lagc = [np.sqrt(a * b) for a, b in AV.LAGS]
    for lab, sel, col in (("elongated", lambda a: a >= 1.25, C2), ("circle-like", lambda a: a < 1.25, C1)):
        obs, exp = AV.clustering([r for r in res if sel(float(mm[r["rel"]]["aspect"]))]); ax[2].plot(lagc, obs / exp, "o-", color=col, ms=4, label=lab)
    ax[2].axhline(1, color="0.6", lw=0.8, ls=":"); ax[2].set_xscale("log"); ax[2].set_xlabel("time after an escape (t.u.)"); ax[2].set_ylabel("later escapes / independent")
    ax[2].set_title("(c) no cascade in elongated cavities", fontsize=9); ax[2].legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(F, "fig_shapes_mech.png"), dpi=200); plt.close(fig)
if __name__ == "__main__":
    lock(); shapes(); ejection()
