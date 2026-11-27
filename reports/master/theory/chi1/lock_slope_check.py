#!/usr/bin/env python3
"""Direct check of the lock slope (user, 2026-11-27): tau_j = s phi_j + c with s = 2k - 2?
Per frame (second half of each run), wall +1/2 (within 4 cells of the wall): azimuth phi_j and tilt tau_j = tail direction - outward normal (= phi_j on a
circle). For a trial slope s, the frame's lock coherence is L(s) = |<exp(i (tau_j - s phi_j))>_j| (1 = every defect obeys tau = s phi + c with one c).
Averaged over frames; scanned s = 0 .. 3k in steps of 0.02. Only integer s gives a single-valued relation on the circle; the scan shows where the peak is.
Data: active MPCD circles R 31.83 (chi 0.5, lambda 2, zeta 0.07), all k, all seeds. Also: an unwrapped circular-linear fit of tau_j vs phi_j per frame.
-> lock_slope_check.txt, fig_lock_slope_check.pdf/.png"""
import glob, json, os, sys, warnings, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa
S_STEP = 0.02


def run_one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
        t, fr = SB.read_with_angle(run, Q, 0.5); ss = np.arange(0, 3 * k + 1e-9, S_STEP); acc = np.zeros_like(ss); n = 0; resid = []
        for f in fr:
            if not len(f): continue
            zs = f[:, 0] + 1j * f[:, 1]; w = (f[:, 2] > 0) & (R - np.abs(zs) < 4)
            if w.sum() < 3: continue
            phi = np.angle(zs[w]); tau = np.angle(np.exp(1j * (f[w, 3] - phi)))
            acc += np.abs(np.exp(1j * (tau[None, :] - ss[:, None] * phi[None, :])).mean(1)); n += 1
            c = np.angle(np.mean(np.exp(1j * (tau - (2 * k - 2) * phi)))); resid += list(np.angle(np.exp(1j * (tau - (2 * k - 2) * phi - c))))
        if n < 10: return None
        return dict(k=k, ss=ss, L=acc / n, n=n, resid=np.degrees(np.array(resid)))
    except Exception as ex:
        return dict(err=repr(ex)[:150], run=run)


if __name__ == "__main__":
    runs = [os.path.dirname(p) for p in sorted(glob.glob(SB.ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/defects.dat.gz"))
            if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(run_one, runs, chunksize=1)
    ok = [r for r in res if r and "err" not in r]; ks = sorted({r["k"] for r in ok})
    L = ["LOCK SLOPE CHECK (active MPCD circles R 31.83): L(s) = <|<exp(i (tau_j - s phi_j))>_j|>_frames",
         "k | runs | frames | best s (scan) | 2k-2 | L(k-1) | L(2k-2) | L(2k) | L(2k-1) | residual about s = 2k-2: circular sd (deg)"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8)); cmap = plt.get_cmap("viridis")
    for i, k in enumerate(ks):
        g = [r for r in ok if r["k"] == k]; ss = g[0]["ss"]; Lm = np.mean([r["L"] for r in g], 0); nf = sum(r["n"] for r in g)
        at = lambda s: Lm[int(round(s / S_STEP))]; best = ss[np.argmax(Lm)]; res_ = np.concatenate([r["resid"] for r in g])
        csd = np.degrees(np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * np.radians(res_)))))))
        L.append(f"{k:2d} | {len(g):2d} | {nf:6d} | {best:5.2f} | {2 * k - 2:2d} | {at(k - 1):.3f} | {at(2 * k - 2):.3f} | {at(2 * k):.3f} | {at(2 * k - 1):.3f} | {csd:5.1f}")
        c = cmap(i / max(len(ks) - 1, 1)); axs[0].plot(ss / (2 * k - 2), Lm, color=c, lw=1.5, label=f"k = {k}")
    for e in [r for r in res if r and "err" in r][:2]: L.append(f"  ERROR {e['run'][-40:]}: {e['err']}")
    axs[0].axvline(1, color="k", ls="--", lw=1); axs[0].axvline(0.5, color="0.6", ls=":", lw=1)
    axs[0].set_xlabel("Trial slope $s$ / $(2k-2)$"); axs[0].set_ylabel("Lock coherence $L(s)$"); axs[0].set_xlim(0, 1.6)
    axs[0].set_title("Which slope makes $\\tau_j - s\\,\\phi_j$ the same for all wall defects?", fontsize=11.5)
    axs[0].text(1.02, 0.95, "$s = 2k-2$", fontsize=10.5); axs[0].text(0.52, 0.95, "$s = k-1$", fontsize=10.5, color="0.45"); axs[0].legend(fontsize=8, ncol=2, frameon=False, loc="center right")
    k0 = 6; g = [r for r in ok if r["k"] == k0]; res_ = np.concatenate([r["resid"] for r in g])
    axs[1].hist(res_, bins=72, range=(-180, 180), color="#1f5fa8"); axs[1].set_xlabel("$\\tau_j - (2k-2)\\phi_j - c$ (deg), k = 6"); axs[1].set_ylabel("Wall-defect frames")
    axs[1].set_title("Residual about the lock, $s = 2k-2$", fontsize=11.5)
    fig.text(0.5, -0.04, "Per frame: azimuth $\\phi_j$ of each wall +½ and its tilt $\\tau_j$ (tail direction minus outward normal); $c$ fitted per frame. Active MPCD circles, R = 31.8, all seeds.",
             ha="center", fontsize=10, color="0.3")
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_lock_slope_check.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_lock_slope_check.png"), dpi=170, bbox_inches="tight")
    open(os.path.join(H, "lock_slope_check.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
