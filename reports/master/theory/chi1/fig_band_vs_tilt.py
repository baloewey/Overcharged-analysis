#!/usr/bin/env python3
"""Layer tilt vs twice the band angle, as a function of k, for six setups (user, 2026-11-24). Band = ring-averaged director angle from the local flow
direction over depths 2-3 cells (complex mean, signed with the circulation); tilt = coherent layer tilt. Circles R 31.83 (MPCD chi 1: R 21.6 = same
R/l_a 1.43). Mean +- s.e. over seeds.
  MPCD chi 0.5: band from curvature/curvature_profiles.npz (1-cell window starting at y = 2), tilt from chi1_compare.one (frames with C > 0.9, 2nd half).
  MPCD chi 1:   band from chi1/chi1_band.npz (same window), tilt from chi1_compare.one.
  Continuum (runs/continuum_chi, laminar zeta 1 / turbulent zeta 8, chi 1 / 0.5): band from the summed profiles (bins [2, 2.5) + [2.5, 3)) over frames with
  |circ| > 0.3 and t >= 500; tilt = mean |tilt| over those frames with coherence > 0.5.
-> fig_band_vs_tilt.pdf/.png, fig_band_vs_tilt.txt"""
import glob, json, os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); import chi1_compare as CC  # noqa
RUNS = "/home/baloewe/Physics/Overcharged/runs"; TH = os.path.join(H, "..")


def mpcd_sets():
    out = {}
    Z = np.load(os.path.join(TH, "curvature", "curvature_profiles.npz")); i2 = int(np.argmin(np.abs(Z["YC"] - 2.25)))
    sel = np.abs(Z["R"] - 31.83) < 0.5; names = [RUNS + "/cavity/server_main/main_act0.07/" + r for r in Z["runs"][sel]]
    with Pool(12) as p: T = p.map(CC.one, names)
    out["MPCD, chi = 0.5"] = [(int(k), 2 * a, t["tilt"]) for k, a, t in zip(Z["k"][sel], Z["alpha"][sel][:, i2], T)]
    B = np.load(os.path.join(H, "chi1_band.npz")); j2 = int(np.argmin(np.abs(B["YC"] - 2.25)))
    names = [RUNS + "/" + r if not r.startswith("/") else r for r in B["runs"]]
    names = [n if os.path.isdir(n) else RUNS + "/chi1_matched/" + os.path.basename(n) for n in names]
    with Pool(12) as p: T = p.map(CC.one, names)
    out["MPCD, chi = 1"] = [(int(k), 2 * a, t["tilt"]) for k, a, t in zip(B["k"], B["alpha"][:, j2], T)]
    return out


def cont_sets():
    out = {}
    for f in sorted(glob.glob(RUNS + "/continuum_chi/z*_k*_chi*_s*.json")):
        J = json.load(open(f)); a = J["args"]; ob = [o for o in J["obs"] if o["t"] >= 500 and abs(o.get("circ", 0)) > 0.3 and "pZr" in o]
        if not ob: continue
        Zs = np.sum([np.array(o["pZr"][4:6]) + 1j * np.array(o["pZi"][4:6]) for o in ob]); band = np.degrees(0.5 * np.angle(Zs))
        tl = [abs(o["tilt"]) for o in ob if "tilt" in o and o.get("coh", 0) > 0.5]
        lab = f"Continuum {'laminar' if a['zeta'] < 3 else 'turbulent'}, chi = {a['chi']:g}"
        out.setdefault(lab, []).append((a["k"], 2 * band, np.mean(tl) if tl else np.nan))
    return out


if __name__ == "__main__":
    D = {**cont_sets(), **mpcd_sets()}
    order = ["Continuum laminar, chi = 1", "Continuum laminar, chi = 0.5", "Continuum turbulent, chi = 1", "Continuum turbulent, chi = 0.5", "MPCD, chi = 1", "MPCD, chi = 0.5"]
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.2), sharey=True); L = ["k | 2 x band (mean +- s.e.) | tilt (mean +- s.e.) | seeds"]
    for ax, lab, let in zip(axs.ravel(), order, "abcdef"):
        rows = D.get(lab, []); ks = sorted({r[0] for r in rows}); L.append(f"\n{lab}")
        mb, eb, mt, et = [], [], [], []
        for k in ks:
            b = np.array([r[1] for r in rows if r[0] == k]); t = np.array([r[2] for r in rows if r[0] == k]); t = t[np.isfinite(t)]
            se = lambda v: v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
            mb.append(b.mean()); eb.append(se(b)); mt.append(t.mean() if len(t) else np.nan); et.append(se(t) if len(t) else 0.0)
            L.append(f"  k {k:2d} | {b.mean():5.1f} +- {se(b):.1f} | {mt[-1]:5.1f} +- {et[-1]:.1f} | {len(b)}")
        ax.axhline(60, color="0.6", ls="--", lw=1.2, label=r"$2\theta_L = \arccos(1/\lambda)$")
        ax.errorbar(ks, mt, yerr=et, fmt="o-", color="#c0392b", ms=8, lw=2, capsize=4, label="Layer tilt")
        ax.errorbar(ks, mb, yerr=eb, fmt="s--", color="#1f5fa8", mfc="white", mew=2, ms=8, lw=2, capsize=4, label=r"$2\times$ band angle")
        ax.set_title(f"({let}) {lab.replace('chi', chr(967))}", fontsize=13); ax.set_xlabel("Winding k"); ax.set_ylim(25, 75)
        if ks: ax.set_xlim(min(ks) - 1, max(ks) + 1)
    for ax in axs[:, 0]: ax.set_ylabel("Angle (deg)")
    axs[0, 0].legend(loc="lower left", fontsize=10, frameon=False)
    fig.text(0.5, -0.01, r"Band: director angle from the flow direction, ring-averaged over depths 2–3 cells. Circles at $R/\ell_a$ = 1.43; mean $\pm$ s.e. over seeds.",
             ha="center", fontsize=11, color="0.3")
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_band_vs_tilt.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_band_vs_tilt.png"), dpi=200, bbox_inches="tight")
    open(os.path.join(H, "fig_band_vs_tilt.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
