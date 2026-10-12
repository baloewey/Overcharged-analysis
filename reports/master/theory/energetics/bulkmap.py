#!/usr/bin/env python3
"""Where does the bulk charge sit in elongated cavities? (user, 2026-10-12). Shape sweep, second half, every 5th frame. Defects deeper than 4 cells, normalised coordinates
u = x/A, v = y/B, elliptic radius rho = sqrt(u^2 + v^2). Density per unit normalised area of +1/2, -1/2 and net charge in rho bins and in |u| bins (centre vs tips), by aspect.
-> bulkmap.txt, fig_bulkmap.png"""
import os, sys, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
RB = np.linspace(0, 1, 11); UB = np.linspace(0, 1, 11)
def one(run):
    try:
        J = __import__("json").load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; ainv = bc.get("aInv", [1, 1, 1]); B = float(bc["R"]); A = B * float(ainv[0] or 1)
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5); out = []; nf = 0
        for f in fr[len(fr) // 2::5]:
            nf += 1
            if not len(f): continue
            d = tree.query(f[:, :2])[0]; m = d > 4
            for x, y, q in f[m, :3]: out.append((x / A, y / B, np.sign(q)))
        return (A / B, np.sqrt(A * B) / 22.21, nf, np.array(out).reshape(-1, 3))
    except Exception as ex:
        return (None, run, str(ex))
if __name__ == "__main__":
    SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; runs = sorted(glob.glob(SW + "/ellipse/*/*/*/s0*") + glob.glob(SW + "/circle/circle/*/*/s0*"))
    with Pool(int(os.environ.get("NPROC", 4))) as pool: res = [r for r in pool.map(one, runs, chunksize=2) if r[0] is not None]
    L = ["Bulk defects (> 4 cells from the wall) per frame per unit normalised area, by elliptic radius rho (0 centre, 1 wall) and by |x|/A; net = (+1/2 count - -1/2 count)/2"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4)); cm = plt.get_cmap("viridis"); ASP = sorted(set(round(r[0], 2) for r in res))
    for ia, a in enumerate(ASP):
        S = [r for r in res if round(r[0], 2) == a]; nf = sum(r[2] for r in S); X = np.concatenate([r[3] for r in S if len(r[3])])
        rho = np.hypot(X[:, 0], X[:, 1]); au = np.abs(X[:, 0]); sp, sm = X[:, 2] > 0, X[:, 2] < 0
        ar = np.pi * np.diff(RB ** 2)                     # normalised area of each rho annulus
        npos = np.histogram(rho[sp], RB)[0] / nf / ar; nneg = np.histogram(rho[sm], RB)[0] / nf / ar; net = (npos - nneg) / 2
        # |u| strips: normalised area of the strip |u| in [u1,u2] inside the unit disk = 2 * integral 2 sqrt(1-u^2) du
        fa = lambda u: u * np.sqrt(1 - u ** 2) + np.arcsin(u); sa = 2 * (fa(UB[1:]) - fa(UB[:-1]))
        upos = np.histogram(au[sp], UB)[0] / nf / sa; uneg = np.histogram(au[sm], UB)[0] / nf / sa; unet = (upos - uneg) / 2
        L.append(f"\naspect {a:.2f} ({len(S)} runs; bulk net charge per frame {(sp.sum() - sm.sum()) / 2 / nf:+.2f}; mean rho of +1/2 {rho[sp].mean():.2f}, uniform 0.67*; mean |x|/A {au[sp].mean():.2f})")
        L.append("  rho bins:  " + " ".join(f"{x:.2f}" for x in RB[1:]) + "\n  net/area:  " + " ".join(f"{x:+.2f}" for x in net) + "\n  +1/2/area: " + " ".join(f"{x:.2f}" for x in npos))
        L.append("  |x|/A bins: " + " ".join(f"{x:.1f}" for x in UB[1:]) + "\n  net/area:   " + " ".join(f"{x:+.2f}" for x in unet))
        ax[0].plot(0.5 * (RB[1:] + RB[:-1]), net, "-o", ms=3, color=cm(ia / max(len(ASP) - 1, 1)), label=f"aspect {a:g}" if a != 1 else "circle")
        ax[1].plot(0.5 * (UB[1:] + UB[:-1]), unet, "-o", ms=3, color=cm(ia / max(len(ASP) - 1, 1)))
    L.append("\n* uniform density inside rho < ~0.8 (depth > 4 cells) gives mean rho of about 0.53 for the full disk; compare between aspects, not to this number.")
    ax[0].set_xlabel("elliptic radius ρ (0 centre, 1 wall)"); ax[0].set_ylabel("bulk net charge per frame per unit normalised area"); ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_xlabel("|x| / A (0 centre, 1 tip)"); ax[1].set_ylabel("bulk net charge per frame per unit normalised area")
    for a_ in ax: a_.axhline(0, color="0.8", lw=0.6); [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.suptitle("Where the bulk charge sits (shape sweep, all sizes, defects > 4 cells from the wall)", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(H, "fig_bulkmap.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "bulkmap.txt"), "w").write("\n".join(L) + "\n")
