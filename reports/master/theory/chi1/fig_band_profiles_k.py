#!/usr/bin/env python3
"""Band angle vs distance to the wall for each winding k, MPCD chi 0.5, circles R 31.83 (user, 2026-11-24). Per run: coherent frames (C > 0.9) with exactly
2k-2 wall +1/2; per 1-cell depth bin [y, y+1), y = 0..11: static-subtracted ring average alpha_corr = alpha_meas - alpha_static + tau_static/2 (angle of the
director from the local flow direction, signed with the layer tilt; band_windows/static_band estimator). Medians over frames, mean over seeds.
Markers: each k's near-field edge d + a/2pi (measured per frame). -> fig_band_profiles_k.pdf/.png, fig_band_profiles_k.txt"""
import glob, gzip, json, os, sys, warnings, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa
YB = np.arange(0, 12, 1.0)


def run_one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = 2 * k - 2
        with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; t, fr = SB.read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        ph = np.arctan2(Y, X); u = R - np.hypot(X, Y); prof, edges, ds = [], [], []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]; zs = f[:, 0] + 1j * f[:, 1]; qs = f[:, 2]; r = np.abs(zs); w = (qs > 0) & (R - r < 4)
            if w.sum() != N or abs(qs.sum() - k) > 1e-6: continue
            ez = np.mean(np.exp(1j * (f[w, 3] - np.angle(zs[w]))))
            if abs(ez) < 0.9: continue
            sg = np.sign(np.angle(ez)) or 1.0; d = np.median(R - r[w]); edges.append(d + (R - d) / N); ds.append(d)
            dd = D[D[:, 0] == tt]; th2 = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th2[ix, iy] = 2 * np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]; th2[u <= 0] = np.nan
            th2s, _, _ = SB.TS.fields(X, Y, list(zs), list(qs), R)
            nt = SB.TS.natural_tilts(list(zs), list(qs), R, len(zs)); ts_ = sg * np.angle(np.mean(np.exp(1j * nt[np.flatnonzero(w)])))
            row = []
            for y in YB:
                b = (u >= y) & (u < y + 1) & np.isfinite(th2)
                am = np.angle(np.sum(S[b] * np.exp(1j * sg * (th2[b] - 2 * ph[b] - np.pi)))) / 2
                ast = np.angle(np.sum(np.exp(1j * sg * (th2s[b] - 2 * ph[b] - np.pi)))) / 2
                row.append(np.angle(np.exp(1j * 2 * (am - ast + ts_ / 2))) / 2)
            prof.append(row)
        if len(prof) < 3: return None
        return dict(k=k, prof=np.degrees(np.median(np.array(prof), 0)), edge=np.median(edges), d=np.median(ds), n=len(prof))
    except Exception as ex:
        return dict(err=repr(ex)[:120], run=run)


if __name__ == "__main__":
    runs = [os.path.dirname(p) for p in sorted(glob.glob(SB.ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
            if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(run_one, runs, chunksize=1)
    ok = [r for r in res if r and "err" not in r]; ks = sorted({r["k"] for r in ok})
    cmap = plt.get_cmap("plasma"); fig, ax = plt.subplots(figsize=(10, 6.4)); L = ["k | runs | near-field edge | alpha(y) for y = " + " ".join(f"{y + 0.5:g}" for y in YB)]
    ax.axhline(30, color="0.5", ls="--", lw=1.3); ax.text(11.6, 30.6, r"$\theta_L$ = 30°", ha="right", va="bottom", fontsize=11, color="0.4")
    for i, k in enumerate(ks):
        g = [r for r in ok if r["k"] == k]; P = np.mean([r["prof"] for r in g], 0); e = np.mean([r["edge"] for r in g]); c = cmap(0.05 + 0.85 * i / max(len(ks) - 1, 1))
        y = YB + 0.5; dd = np.mean([r["d"] for r in g])
        inside = y < dd + 0.5; near = (~inside) & (y < e); beyond = y >= e - 0.5         # inside the defect row: the subtraction returns ~tau/2 (masked)
        dot = (~inside) & (y < e + 0.5); ax.plot(y[dot], P[dot], ":", color=c, lw=1.6)
        sol = (~inside) & beyond; ax.plot(y[sol], P[sol], "-", color=c, lw=2.4, label=f"k = {k}")
        ax.plot(e, np.interp(e, y, P), "o", color=c, ms=10, mec="k", mew=1.2, zorder=5)
        L.append(f"{k:2d} | {len(g):2d} | {e:5.2f} | " + " ".join(f"{v:5.1f}" for v in P))
    ax.plot([], [], "o", color="w", mec="k", mew=1.2, ms=10, label="near-field edge d + a/2π"); ax.plot([], [], ":", color="0.4", lw=1.6, label="inside the near field (between the defect row and the edge)")
    ax.set_xlabel("Distance to the wall (cells)"); ax.set_ylabel("Band angle: director angle from the flow (deg)"); ax.set_xlim(0, 12); ax.set_ylim(-5, 40)
    ax.set_title(r"MPCD, $\chi$ = 0.5, circles $R$ = 31.8 ($R/\ell_a$ = 1.43): band profile for each winding", fontsize=13)
    ax.legend(fontsize=9.5, frameon=False, ncol=2, loc="lower left")
    fig.text(0.5, -0.04, "Static-subtracted ring average in 1-cell depth bins; coherent frames with exactly 2k−2 wall defects; mean over seeds.\nDepths at or inside the defect row (y < d + 0.5) are omitted: there the subtraction returns half the layer tilt, not a flow-aligned band.", ha="center", fontsize=10.5, color="0.3")
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_band_profiles_k.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_band_profiles_k.png"), dpi=200, bbox_inches="tight")
    open(os.path.join(H, "fig_band_profiles_k.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
