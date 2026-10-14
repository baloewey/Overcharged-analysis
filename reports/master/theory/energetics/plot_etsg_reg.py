#!/usr/bin/env python3
"""Shape-sweep energetics split by wall region (etsg_shapes_reg.npy): depleting escapers vs resident controls (matched in |x|/A) at the tips / mid / flanks.
Pooled over aspect >= 1.4 (all sizes); aspect 1.2 shown separately (tips dominate there). -> fig_etsg_reg.png, etsg_reg.txt"""
import os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); a = np.load(os.path.join(H, "etsg_shapes_reg.npy"), allow_pickle=True)
kind = a[:, 1].astype(str); lag = a[:, 2].astype(int); F0, E6, ct, asp = (a[:, i].astype(float) for i in (3, 5, 6, 7)); rg = a[:, 9].astype(int)
RN = ("tip", "mid", "flank"); COL = ("#c2452d", "#7a3b8f", "#1f5fa8")
def ms(x): x = x[np.isfinite(x)]; return (x.mean(), x.std() / np.sqrt(len(x)), len(x)) if len(x) else (np.nan, np.nan, 0)
fig, ax = plt.subplots(1, 3, figsize=(14, 4.6)); L = []
for sel_name, sel in (("aspect >= 1.4", asp >= 1.35), ("aspect 1.2", np.abs(asp - 1.2) < 0.05)):
    for r in range(3):
        for q, (Y, nm) in enumerate(((ct, "cos tau"), (F0, "F0"), (E6, "E6"))):
            for kd, ls in (("esc", "-"), ("ctl", "--")):
                m0 = sel & (kind == kd) & (rg == r); ls_ = sorted(set(lag[m0]))
                if not ls_: continue
                st = [ms(Y[m0 & (lag == l)]) for l in ls_]
                L.append(f"{sel_name:13s} {RN[r]:5s} {kd} {nm:7s}: " + "  ".join(f"{2*l:+d}:{x[0]:+.2f}±{x[1]:.2f}(n{x[2]})" for l, x in zip(ls_, st)))
                if sel_name == "aspect >= 1.4":
                    ax[q].errorbar(2 * np.array(ls_), [x[0] for x in st], [x[1] for x in st], color=COL[r], ls=ls, marker="o" if kd == "esc" else "s", ms=3, capsize=2, lw=1.6 if kd == "esc" else 1,
                                   label=f"{RN[r]} ({st[-1][2] if kd == 'esc' else ''} escapes)" if (kd == "esc" and q == 0) else None)
for q, yl in enumerate(("cos τ", "F₀ (kBT/cell)", "E₆ (kBT)")):
    ax[q].set_xlabel("time relative to departure (t.u.)"); ax[q].set_ylabel(yl)
    ax[q].set_title(("(a) tilt projection (1 = head into the bulk)", "(b) elastic force along the inward normal", "(c) elastic cost of moving 6 cells in")[q], fontsize=10); ax[q].axvline(0, color="0.6", lw=0.8, ls=":")
    for s_ in ("top", "right"): ax[q].spines[s_].set_visible(False)
fig.legend(*ax[0].get_legend_handles_labels(), fontsize=8, frameon=False, ncol=3, loc="lower center", title="Solid: escaper, dashed: control", title_fontsize=8)
fig.tight_layout(rect=(0, 0.12, 1, 1)); fig.savefig(os.path.join(H, "fig_etsg_reg.png"), dpi=200)
open(os.path.join(H, "etsg_reg.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
