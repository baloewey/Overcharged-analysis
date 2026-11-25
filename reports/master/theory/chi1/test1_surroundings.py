#!/usr/bin/env python3
"""Test 1 of the working picture (user, 2026-11-25): what sets the field just outside each wall defect's core?
For a +1/2, theta ~ (1/2) arg(z - z_j) + c_j near the core and the tail lies at 2 c_j: "tilt = 2 x surroundings" is exact AT the core. So we measure how the
background phase c_j(r) varies with distance r from the core: per wall defect (coherent frames, exactly 2k-2 wall +1/2), the field is rotated into the defect's
local frame (outward normal along x), the defect's OWN singular part (1/2) arg(z - z_j) is removed (neighbours kept), and c_j(r) = (1/2) arg < S exp(2i(theta' -
(1/2) arg(z' - z_j'))) > over cells at distance [r, r+1) from the core that lie DEEPER than the defect (bulk side; the wall is never sampled). 2 c_j(r) is the tilt
the defect would have if locked to its surroundings at distance r; signed with the layer tilt. Prediction (working picture): crowded -> 2 x band quickly;
sparse -> stays near 2 theta_L (near-wall simple shear) out to several cells.
Sets: MPCD chi 0.5 (R 31.8) and chi 1 (R 21.6). Also the layer tilt and 2 x band just outside the near field (band_windows.npy) as references.
-> test1_surroundings.npy, test1_surroundings.txt, fig_test1_surroundings.pdf/.png"""
import glob, gzip, json, os, sys, warnings, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa
RS = np.arange(1.0, 8.0, 1.0)


def run_one(args):
    run, stride = args
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = 2 * k - 2
        with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])][::stride]; t, fr = SB.read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        Z = (X + 1j * Y).ravel(); RR = np.abs(Z); prof, taus = [], []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]; zs = f[:, 0] + 1j * f[:, 1]; qs = f[:, 2]; r = np.abs(zs); w = (qs > 0) & (R - r < 4)
            if w.sum() != N or abs(qs.sum() - k) > 1e-6: continue
            ez = np.mean(np.exp(1j * (f[w, 3] - np.angle(zs[w]))))
            if abs(ez) < 0.9: continue
            sg = np.sign(np.angle(ez)) or 1.0
            dd = D[D[:, 0] == tt]; th = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th[ix, iy] = np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]; th = th.ravel(); S = S.ravel(); ok = np.isfinite(th) & (RR < R - 0.5)
            for j in np.flatnonzero(w):
                zj = zs[j]; phj = np.angle(zj); dz = Z - zj; dist = np.abs(dz)
                deeper = RR < np.abs(zj)                                                   # bulk side of the defect
                thp = th - phj; argp = np.angle(dz) - phj                                  # local frame: rotate field and positions by -phi_j
                bg = 2 * thp - argp                                                        # 2 x background phase (own +1/2 removed), headless
                row = []
                for r0 in RS:
                    sel = ok & deeper & (dist >= r0) & (dist < r0 + 1)
                    row.append(np.angle(np.sum(S[sel] * np.exp(1j * bg[sel]))) if sel.sum() >= 4 else np.nan)
                prof.append(sg * np.array(row)); taus.append(sg * np.angle(np.exp(1j * (f[j, 3] - phj))))
        if len(prof) < 20: return None
        P = np.array(prof); return dict(k=k, R=R, prof=np.degrees(np.angle(np.nanmean(np.exp(1j * P), 0))), tau=np.degrees(np.angle(np.mean(np.exp(1j * np.array(taus))))), n=len(prof))
    except Exception as ex:
        return dict(err=repr(ex)[:150], run=run)


if __name__ == "__main__":
    c05 = [(os.path.dirname(p), 1) for p in sorted(glob.glob(SB.ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
           if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    c1 = [(os.path.dirname(p), 5) for p in sorted(glob.glob(SB.ROOT + "/runs/chi1_matched/R1.43_k*_s*/directorfield.dat.gz"))]
    with Pool(int(os.environ.get("NPROC", 12))) as p: r05 = p.map(run_one, c05, chunksize=1); r1 = p.map(run_one, c1, chunksize=1)
    BW = np.load(os.path.join(H, "band_windows.npy"), allow_pickle=True).item()
    L = ["TEST 1: background phase around each wall defect vs distance from the core (2 c(r), deg; own +1/2 removed, bulk side only, local frame)",
         "set | k | defect-frames | defect tilt | 2c(r) for r = " + " ".join(f"{r:g}-{r + 1:g}" for r in RS) + " | 2 x band outside near field | 2 theta_L = 60"]
    out = {}; fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
    for ax, (lab, rr) in zip(axs, (("MPCD, chi = 0.5", r05), ("MPCD, chi = 1", r1))):
        ok = [r for r in rr if r and "err" not in r]; ks = sorted({r["k"] for r in ok}); cmap = plt.get_cmap("plasma")
        for e in [r for r in rr if r and "err" in r][:2]: L.append(f"  ERROR {e['run'][-40:]}: {e['err']}")
        bands = {row["k"]: row["ec"] for row in BW.get(lab, [])}
        for i, k in enumerate(ks):
            g = [r for r in ok if r["k"] == k]; P = np.degrees(np.angle(np.mean([np.exp(1j * np.radians(r["prof"])) for r in g], 0)))
            tau = np.mean([r["tau"] for r in g]); nn = sum(r["n"] for r in g); c = cmap(0.05 + 0.85 * i / max(len(ks) - 1, 1))
            out.setdefault(lab, []).append(dict(k=k, tau=tau, prof=P, band=bands.get(k, np.nan)))
            ax.plot(np.r_[0, RS + 0.5], np.r_[tau, P], "o-", color=c, lw=2, ms=5, label=f"k = {k}")
            if np.isfinite(bands.get(k, np.nan)): ax.plot([7.4, 7.9], [bands[k]] * 2, "-", color=c, lw=4)
            L.append(f"{lab:16s} | {k:2d} | {nn:6d} | {tau:5.1f} | " + " ".join(f"{v:5.1f}" for v in P) + f" | {bands.get(k, np.nan):5.1f} | 60")
        ax.axhline(60, color="0.5", ls="--", lw=1.3); ax.text(0.1, 61, r"$2\theta_L$", color="0.4", fontsize=11)
        ax.set_title(f"{lab.replace('chi', chr(967))}: surroundings of each wall defect", fontsize=12.5); ax.set_xlabel("Distance from the defect core (cells; 0 = defect tilt)")
        ax.legend(fontsize=8.5, frameon=False, ncol=2, loc="lower left"); ax.set_xlim(-0.2, 8.0)
    axs[0].set_ylabel(r"$2\,c(r)$: tilt if locked to the surroundings at $r$ (deg)"); axs[0].set_ylim(0, 75)
    fig.text(0.5, -0.03, "Background phase of the director around each wall +½ (its own singular part removed, neighbours kept), cells deeper than the defect only. "
             "Bars at right: 2 × band just outside the near field.", ha="center", fontsize=10.5, color="0.3")
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_test1_surroundings.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_test1_surroundings.png"), dpi=180, bbox_inches="tight")
    np.save(os.path.join(H, "test1_surroundings.npy"), out, allow_pickle=True)
    open(os.path.join(H, "test1_surroundings.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
