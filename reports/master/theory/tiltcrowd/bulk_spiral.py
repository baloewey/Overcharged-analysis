#!/usr/bin/env python3
"""Does the bulk spiral follow the wall-layer tilt? (user, 2026-10-15: a smaller tilt changes the nematic spiral and with it the bulk flow.)
Effective boundary: behind a layer tilted by tau the bulk sees a planar condition rotated by ~ tau/2 -> the bulk +1 is a spiral with angle beta ~ tau/2
(director vs radial direction); its active drive of the circulation ~ sin 2 beta.
Runs: main circle sweep at fixed R = 31.83 (1.43 l_a), k = 3-15 (10 seeds) and 16, 20 (3 seeds) (spacing 2 pi R/(2k-2) = 50 ... 5.3 cells); second half.
Per director frame (every 100 t.u.): beta(r) = (1/2) arg < S exp(2i(theta - phi)) > over annuli (2-cell bins), signed by the circulation sense of the flow
frame at the same time (beta > 0 = the director leans with the circulation); bulk spiral = median over 0.3 R < r < 0.75 R; order |<S e^{2i(theta-phi)}>|.
Also: tilt of the coherent wall layer (C > 0.9) and |<v_theta>| (interior). -> bulk_spiral.txt, fig_bulk_spiral.png"""
import os, sys, glob, gzip, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa


def frames(fn, ncol):
    rows = np.loadtxt((l for l in gzip.open(fn, "rt") if len(l.split()) == ncol and l.split()[0][0].isdigit())); return rows


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    D = frames(run + "/directorfield.dat.gz", 8); F = frames(run + "/flowfield.dat.gz", 7)
    td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]
    edges = np.arange(0, R, 2.0); rc = 0.5 * (edges[1:] + edges[:-1]); prof = []; order = []; vth = []; zs = []
    for tt in td:
        d = D[D[:, 0] == tt]; f = F[np.abs(F[:, 0] - tt) < 1e-6]
        if not len(f): continue
        x = d[:, 1] + 0.5 - Q[0]; y = d[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); phi = np.arctan2(y, x); m = (r < R - 0.5) & (d[:, 7] > 0)
        xf = f[:, 1] + 0.5 - Q[0]; yf = f[:, 2] + 0.5 - Q[1]; rf = np.hypot(xf, yf); mf = rf < R - 1
        vt = (xf * f[:, 5] - yf * f[:, 4]) / np.maximum(rf, 1e-9); circ = np.sign(vt[mf].mean())
        if circ == 0: continue
        vth.append(abs(vt[mf & (rf > 0.3 * R)].mean()))
        z = d[m, 7] * np.exp(2j * (np.arctan2(d[m, 5], d[m, 4]) - phi[m]))
        num, _ = np.histogram(r[m], edges, weights=z.real); numi, _ = np.histogram(r[m], edges, weights=z.imag); cnt, _ = np.histogram(r[m], edges)
        zz = (num + 1j * numi) / np.maximum(cnt, 1)
        prof.append(circ * 0.5 * np.angle(zz)); order.append(np.abs(zz)); zs.append(np.abs(zz) * np.exp(1j * circ * np.angle(zz)))
    P = np.array(prof); O = np.array(order); band = (rc > 0.3 * R) & (rc < 0.75 * R)
    t, fr = read_with_angle(run, Q, 0.5); tl = []
    for fz in fr:
        if not len(fz): continue
        rr = np.hypot(fz[:, 0], fz[:, 1]); w = (fz[:, 2] > 0) & (R - rr < 4)
        if w.sum() < 3: continue
        zt = np.mean(np.exp(1j * (fz[w, 3] - np.arctan2(fz[w, 1], fz[w, 0]))))
        if abs(zt) > 0.9: tl.append(abs(np.angle(zt)))
    return dict(run=run, k=k, R=R, a=2 * np.pi * R / (2 * k - 2), beta=np.degrees(np.median(P[:, band])), beta_prof=np.degrees(np.median(P, 0)), rc=rc / R,
                order=np.median(O[:, band]), vth=np.mean(vth), zprof=np.mean(zs, 0), tilt=np.degrees(np.median(tl)) if len(tl) > 20 else np.nan)


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz")
                  if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.05)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs)
    ks = sorted(set(x["k"] for x in res))
    L = ["BULK SPIRAL vs WALL TILT (main circle sweep, R = 31.83 = 1.43 l_a, second half). beta = director angle from radial in 0.3R < r < 0.75R, signed by circulation (deg).",
         "k | seeds | spacing | wall tilt tau (C>0.9) | tau/2 | bulk beta (mean +- s.e.) | spiral order | sin 2 beta | sin(tau) | |v_theta| interior"]
    G = {}
    for kk in ks:
        S = [x for x in res if x["k"] == kk]; b = np.array([x["beta"] for x in S]); ta = np.nanmean([x["tilt"] for x in S]); v = np.mean([x["vth"] for x in S])
        G[kk] = (S[0]["a"], ta, b.mean(), b.std(ddof=1) / np.sqrt(len(b)) if len(b) > 1 else 0, np.mean([x["order"] for x in S]), v, np.median([x["beta_prof"] for x in S], 0), S[0]["rc"])
        L.append(f"{kk:2d} | {len(S):2d} | {S[0]['a']:5.1f} | {ta:5.1f} | {ta/2:5.1f} | {b.mean():5.1f} +- {G[kk][3]:.1f} | {G[kk][4]:.2f} | {np.sin(np.radians(2*b.mean())):.2f} | {np.sin(np.radians(ta)):.2f} | {v:.3f}")
    L.append("Time- and seed-averaged spiral (circulation-signed complex order, beta in [0, 180) deg; 90 = tangential): beta (order) at r/R = 0.4 / 0.5 / 0.6 / 0.75 / 0.85")
    for kk in ks:
        S = [x for x in res if x["k"] == kk]; Z = np.mean([x["zprof"] for x in S], 0); rc_ = S[0]["rc"]; cells = []
        for rr in (0.4, 0.5, 0.6, 0.75, 0.85):
            i = int(np.argmin(np.abs(rc_ - rr))); cells.append(f"{np.degrees(0.5 * np.angle(Z[i])) % 180:5.1f} ({abs(Z[i]):.2f})")
        L.append(f"  k {kk:2d} (spacing {S[0]['a']:4.1f}, tau/2 {np.nanmean([x['tilt'] for x in S]) / 2:4.1f}): " + "  ".join(cells))
    a = np.array([G[kk][0] for kk in ks]); ta = np.array([G[kk][1] for kk in ks]); be = np.array([G[kk][2] for kk in ks]); vv = np.array([G[kk][5] for kk in ks])
    L.append(f"corr(beta, tau/2) {np.corrcoef(be, ta/2)[0,1]:+.2f}; slope beta vs tau/2 {np.polyfit(ta/2, be, 1)[0]:.2f}; corr(|v_theta|, sin 2 beta) {np.corrcoef(vv, np.sin(np.radians(2*be)))[0,1]:+.2f}")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2)); cols = plt.cm.viridis(np.linspace(0, 0.9, len(ks)))
    ax[0].plot(a, ta / 2, "o-", color="#c2452d", label="half the wall-layer tilt, $\\tau/2$"); ax[0].errorbar(a, be, [G[kk][3] for kk in ks], fmt="s-", color="#1f5fa8", capsize=2, label="bulk spiral angle $\\beta$")
    ax[0].set_xscale("log"); ax[0].set_xticks([5, 7, 10, 15, 20, 30, 50]); ax[0].set_xticklabels(["5", "7", "10", "15", "20", "30", "50"]); ax[0].minorticks_off()
    ax[0].set_xlabel("wall-defect spacing (cells)"); ax[0].set_ylabel("angle (deg)"); ax[0].legend(fontsize=8, frameon=False); ax[0].set_title("(a) Bulk spiral vs wall tilt (R = 1.43 $\\ell_a$)", fontsize=10)
    for kk, c in zip(ks, cols): ax[1].plot(G[kk][7], G[kk][6], color=c, lw=1.2, label=f"k = {kk}")
    ax[1].set_xlabel("r / R"); ax[1].set_ylabel("spiral angle $\\beta(r)$ (deg)"); ax[1].legend(fontsize=6.5, frameon=False, ncol=2); ax[1].set_title("(b) Spiral profiles", fontsize=10)
    ax[2].plot(np.sin(np.radians(2 * be)), vv, "o", color="0.3"); [ax[2].annotate(f"{kk}", (np.sin(np.radians(2 * b_)), v_), fontsize=7, xytext=(3, 2), textcoords="offset points") for kk, b_, v_ in zip(ks, be, vv)]
    ax[2].set_xlabel("active drive of the spiral, sin 2$\\beta$"); ax[2].set_ylabel("interior $|\\langle v_\\theta\\rangle|$"); ax[2].set_title("(c) Does the spiral set the flow? (labels: k)", fontsize=10)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_bulk_spiral.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "bulk_spiral.txt"), "w").write("\n".join(L) + "\n")
