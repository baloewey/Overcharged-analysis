#!/usr/bin/env python3
"""Talk variants of Fig. 3(b) (user, 2026-11-25), same size and style as the f3b crop:
 f3b_lam2.png    : lambda >= 2 only in the main panel (inset unchanged, all lambda)
 f3b_coh.png     : <|tau|> over coherent frames only (C > 0.9), all lambda; fraction of coherent frames written next to the lambda < 2 points
 f3b_signed.png  : per-frame |signed mean tilt| |<tau>|, all lambda
Data: energetics/tiltlaw.txt (all) and energetics/tiltlaw_coherent.npy (coh, signed). Run: ../paper/prl/figs/run.sh-style env (prlstyle on the path)."""
import os, sys, numpy as np
F = "/home/baloewe/Physics/Overcharged/reports/paper/prl/figs"; sys.path.insert(0, F); os.chdir(F)
import importlib.util
spec = importlib.util.spec_from_file_location("f3", F + "/fig3_tiltlaw.py"); f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)
plt = f3.plt; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
Z = np.load(f3.TH + "/energetics/tiltlaw_coherent.npy", allow_pickle=True).item()
for name, mode, lammin, ylab in (("f3b_lam2", "all", 2.0, None), ("f3b_coh", "coh", 0.0, r"Tilt $\langle|\tau|\rangle$, coherent frames (deg)"),
                                 ("f3b_signed", "signed", 0.0, r"Common tilt $|\langle\tau\rangle|$ (deg)")):
    fig, ax = plt.subplots(figsize=(0.367 * f3.PAGE, 2.25)); fig.subplots_adjust(left=0.2, right=0.97, bottom=0.17, top=0.93)
    f3.tiltlaw(ax, mode=mode, lammin=lammin)
    if ylab: ax.set_ylabel(ylab)
    if mode == "coh":
        for (k, l), d in Z.items():
            if l < 2 and np.isfinite(d["coh"]):
                ax.annotate(f"{100 * d['fcoh']:.0f}%", (l, d["coh"]), xytext=(4, -2 if k != 4 else 4), textcoords="offset points", fontsize=5.2, color=f3.KCOL[k])
        ax.text(0.98, 0.56, "% = share of frames\ncoherent (C > 0.9)", transform=ax.transAxes, fontsize=5.2, ha="right", va="center", color="0.35")
    for ext in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    print("wrote", name)
