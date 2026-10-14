#!/usr/bin/env python3
"""Does the bulk overcharge 'dope' the dynamical defect population? (user, 2026-10-14). Winding channels at fixed geometry (W 40, L 192), imposed charge 2m varied
(channel_qscan m 3, 9, 12 + channel m 0, 6). Hypothesis (law of mass action, as n p = n_i^2 in a semiconductor): pairs are created at a rate per area and annihilate
at a rate proportional to n+ n-, so n+ n- = n0^2 and an excess delta = n+ - n- of +1/2 suppresses -1/2 as n- = (sqrt(delta^2 + 4 n0^2) - delta)/2.
Per run (second half; frames every 2 t.u.), BULK = more than 4 cells from both walls:
  * number densities n+, n- of bulk +1/2 and -1/2 (per cell^2), excess delta = n+ - n-, product <n+ n->; their fluctuations;
  * pair creation and annihilation events (new / vanishing +1/2 and -1/2 within 3 cells; frame-to-frame links < 2.5 cells, periodic in x), rates per bulk area
    per 1000 t.u., and where they happen across the width;
  * like- and unlike-charge pair correlations of bulk defects at short range (count within r < 6 cells relative to a uniform distribution at the same density);
  * mean speed of bulk +1/2 and -1/2 (displacement over one frame).
Fit: one parameter n0 for the doping law over all runs; also n- vs delta per run. Partial runs (no DONE) are flagged.
Usage: doping.py RUNDIR [...] -> doping.txt, fig_doping.png (+pdf with the house style)"""
import sys, os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); BAND = 4.0


def link(f, g, L):
    """nearest same-sign defect of g for each defect of f within 2.5 cells (periodic in x); -1 if none"""
    out = -np.ones(len(f), int)
    for i in range(len(f)):
        s = np.flatnonzero(np.sign(g[:, 2]) == np.sign(f[i, 2]))
        if not len(s): continue
        dx = (g[s, 0] - f[i, 0] + L / 2) % L - L / 2; dd = np.hypot(dx, g[s, 1] - f[i, 1]); k = np.argmin(dd)
        if dd[k] < 2.5: out[i] = s[k]
    return out


def pairs(f, idx, sgn_a, sgn_b, L):
    """+/- pairs among defects idx (unlinked ones) within 3 cells; returns list of mean y"""
    used, ys = set(), []
    for i in idx:
        if np.sign(f[i, 2]) != sgn_a or i in used: continue
        for j in idx:
            if j not in used and np.sign(f[j, 2]) == sgn_b and np.hypot((f[j, 0] - f[i, 0] + L / 2) % L - L / 2, f[j, 1] - f[i, 1]) < 3:
                ys.append(0.5 * (f[i, 1] + f[j, 1])); used |= {i, j}; break
    return ys


def analyse(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; W = y1 - y0
    m = J["BC"][2].get("customAnch", 0); done = os.path.exists(run + "/DONE")
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); dt = float(np.median(np.diff(t))); Ab = L * (W - 2 * BAND)
    npl, nmi, prod, nuc, ann, sp, sm, close = [], [], [], [], [], [], [], np.zeros(3); exp_close = np.zeros(3); prev = None; nfr = 0
    for f in fr:
        if not len(f): prev = None; continue
        nfr += 1; bulk = (f[:, 1] - y0 >= BAND) & (y1 - f[:, 1] >= BAND)
        p = (f[:, 2] > 0) & bulk; q = (f[:, 2] < 0) & bulk; npl.append(p.sum()); nmi.append(q.sum()); prod.append(p.sum() * q.sum())
        # short-range pair correlations among bulk defects: ++, --, +- within r < 6 (periodic in x), vs uniform expectation pi r^2 n / area
        P, M = f[p, :2], f[q, :2]
        def cnt(A, B, same):
            if len(A) == 0 or len(B) == 0: return 0
            dx = (A[:, None, 0] - B[None, :, 0] + L / 2) % L - L / 2; dd = np.hypot(dx, A[:, None, 1] - B[None, :, 1])
            c = (dd < 6).sum() - (len(A) if same else 0); return c / 2 if same else c
        close += [cnt(P, P, True), cnt(M, M, True), cnt(P, M, False)]
        a6 = np.pi * 36 / Ab; exp_close += [len(P) * (len(P) - 1) / 2 * a6, len(M) * (len(M) - 1) / 2 * a6, len(P) * len(M) * a6]
        if prev is not None and len(prev):
            fw = link(prev, f, L); bw = link(f, prev, L)
            for i, j in enumerate(fw):
                if j >= 0 and prev[i, 1] - y0 >= BAND and y1 - prev[i, 1] >= BAND:
                    v = np.hypot((f[j, 0] - prev[i, 0] + L / 2) % L - L / 2, f[j, 1] - prev[i, 1]) / dt
                    (sp if prev[i, 2] > 0 else sm).append(v)
            ann += pairs(prev, [i for i in range(len(prev)) if fw[i] < 0], 1, -1, L)
            nuc += pairs(f, [i for i in range(len(f)) if bw[i] < 0], 1, -1, L)
        prev = f
    T = nfr * dt; npl, nmi, prod = np.array(npl, float), np.array(nmi, float), np.array(prod, float)
    inb = lambda ys: np.array([(y - y0 >= BAND) and (y1 - y >= BAND) for y in ys], bool)
    nb, ab = inb(nuc), inb(ann)
    return dict(run=os.path.relpath(run, "/home/baloewe/Physics/Overcharged/runs"), m=m, W=W, L=L, done=done, T=T, Ab=Ab,
                np=npl.mean() / Ab, nm=nmi.mean() / Ab, prod=prod.mean() / Ab ** 2, np_sd=npl.std() / Ab, nm_sd=nmi.std() / Ab,
                nuc=1000 * nb.sum() / (Ab * T) if len(nuc) else 0.0, ann=1000 * ab.sum() / (Ab * T) if len(ann) else 0.0,
                nuc_bulkfrac=nb.mean() if len(nuc) else np.nan, ann_bulkfrac=ab.mean() if len(ann) else np.nan, n_nuc=len(nuc), n_ann=len(ann),
                g=close / np.maximum(exp_close, 1e-9), vp=np.mean(sp) if sp else np.nan, vm=np.mean(sm) if sm else np.nan)


if __name__ == "__main__":
    runs = [os.path.abspath(a) for a in sys.argv[1:] if os.path.isdir(a)]
    R = [analyse(r) for r in runs]; R.sort(key=lambda r: (r["m"], r["run"]))
    d = np.array([r["np"] - r["nm"] for r in R]); nm = np.array([r["nm"] for r in R]); npl = np.array([r["np"] for r in R])
    # one-parameter fit of the doping law n- = (sqrt(d^2 + 4 n0^2) - d)/2 over runs with defects
    from scipy.optimize import minimize_scalar
    ok = (npl + nm) > 0
    f = lambda n0: np.sum((nm[ok] - (np.sqrt(d[ok] ** 2 + 4 * n0 ** 2) - d[ok]) / 2) ** 2)
    n0 = minimize_scalar(f, bounds=(1e-6, 1e-2), method="bounded").x
    L_ = ["DOPING TEST in winding channels (bulk = more than 4 cells from both walls; second half of each run). Densities per cell^2 (x1e-3), rates per bulk area per 1000 t.u. (x1e-3).",
          f"Mass-action fit n- = (sqrt(delta^2 + 4 n0^2) - delta)/2 over all runs: n0 = {n0 * 1e3:.3f}e-3 per cell^2",
          "run | m (imposed charge 2m) | done | n+ | n- | excess delta | <n+ n-> (x1e-6) | doping-law n- | creation | annihilation | creation in bulk | g++ g-- g+- (r<6) | speed +1/2 -1/2"]
    for r in R:
        dd = r["np"] - r["nm"]; pred = (np.sqrt(dd ** 2 + 4 * n0 ** 2) - dd) / 2
        L_.append(f"{r['run']:28s} | {r['m']:2d} | {'yes' if r['done'] else 'PARTIAL':7s} | {r['np']*1e3:6.3f} | {r['nm']*1e3:6.3f} | {dd*1e3:+6.3f} | {r['prod']*1e6:6.3f} | {pred*1e3:6.3f} | "
                  f"{r['nuc']*1e3:7.3f} | {r['ann']*1e3:7.3f} | {r['nuc_bulkfrac']:.2f} | {r['g'][0]:.2f} {r['g'][1]:.2f} {r['g'][2]:.2f} | {r['vp']:.3f} {r['vm']:.3f}")
    print("\n".join(L_)); open(os.path.join(H, os.environ.get("OUT", "doping.txt")), "w").write("\n".join(L_) + "\n")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4)); cm = {0: "0.5", 3: "#4a86b8", 6: "#7aa36b", 9: "#e08a3c", 12: "#c2452d"}
    for r in R:
        dd = (r["np"] - r["nm"]) * 1e3; c = cm.get(r["m"], "k"); mk = "o" if r["done"] else "^"
        ax[0].plot(dd, r["nm"] * 1e3, mk, color=c, ms=7); ax[0].plot(dd, r["np"] * 1e3, mk, color=c, ms=7, mfc="w")
        ax[1].plot(dd, r["prod"] * 1e6, mk, color=c, ms=7)
        ax[2].plot(dd, r["nuc"] * 1e3, mk, color=c, ms=7); ax[2].plot(dd, r["ann"] * 1e3, "s", color=c, ms=5, mfc="w")
    x = np.linspace(0, max(d.max() * 1e3, 0.1) * 1.05, 100)
    ax[0].plot(x, (np.sqrt(x ** 2 + 4 * (n0 * 1e3) ** 2) - x) / 2, "k--", lw=1, label=f"mass action, $n_0$ = {n0*1e3:.2f}$\\times10^{{-3}}$")
    ax[0].plot(x, (np.sqrt(x ** 2 + 4 * (n0 * 1e3) ** 2) + x) / 2, "k:", lw=1)
    ax[1].axhline((n0 * 1e3) ** 2, color="k", ls="--", lw=1, label="$n_0^2$ (fit)")
    ax[0].set_xlabel("excess density $n_+-n_-$ ($10^{-3}$ per cell$^2$)"); ax[0].set_ylabel("bulk density ($10^{-3}$ per cell$^2$)")
    ax[0].set_title("(a) $-1/2$ (filled) and $+1/2$ (open) densities", fontsize=10); ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_xlabel("excess density $n_+-n_-$ ($10^{-3}$ per cell$^2$)"); ax[1].set_ylabel(r"$\langle n_+n_-\rangle$ ($10^{-6}$ per cell$^4$)")
    ax[1].set_title("(b) mass-action product", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
    ax[2].set_xlabel("excess density $n_+-n_-$ ($10^{-3}$ per cell$^2$)"); ax[2].set_ylabel("events per bulk area per 1000 t.u. ($10^{-3}$)")
    ax[2].set_title("(c) creation (filled), annihilation (open)", fontsize=10)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=cm[k], label=f"m = {k}") for k in sorted(cm)] + [Line2D([], [], marker="^", ls="", color="0.3", label="partial run")],
               fontsize=8, frameon=False, ncol=6, loc="lower center")
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(os.path.join(H, os.environ.get("FIG", "fig_doping.png")), dpi=200)
