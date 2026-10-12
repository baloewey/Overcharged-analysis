#!/usr/bin/env python3
"""Shape-sweep energetics of depleting escapes (etsg_shapes.npy from energetics_general.py): by aspect ratio, escaper (lags -20..+2 frames, 2 t.u. each) vs resident control
(-20, -10, 0): cos tau, elastic restoring force F0, cost to move 6 cells in E6. Means +- s.e. (all sizes pooled; size split printed). -> fig_etsg_shapes.png, etsg_shapes.txt"""
import os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); a = np.load(os.path.join(H, "etsg_shapes.npy"), allow_pickle=True)
kind = a[:, 1].astype(str); lag = a[:, 2].astype(int); F0, E6, ct, asp, req = (a[:, i].astype(float) for i in (3, 5, 6, 7, 8))
asp = np.round(asp, 2); req = np.round(req, 1); ASP = sorted(set(asp)); cmap = plt.get_cmap("viridis")
col = {s: cmap(0.08 + 0.8 * i / (len(ASP) - 1)) for i, s in enumerate(ASP)}
def ms(x): x = x[np.isfinite(x)]; return (x.mean(), x.std() / np.sqrt(len(x)), len(x)) if len(x) else (np.nan, np.nan, 0)
fig, ax = plt.subplots(1, 3, figsize=(13, 4)); L = []
for s in ASP:
    for q, (Y, nm) in enumerate(((ct, "cos tau"), (F0, "F0"), (E6, "E6"))):
        for kd, ls, mk in (("esc", "-", "o"), ("ctl", "--", "s")):
            ls_ = sorted(set(lag[(kind == kd) & (asp == s)])); st = [ms(Y[(kind == kd) & (asp == s) & (lag == l)]) for l in ls_]
            m = np.array([x[0] for x in st]); e = np.array([x[1] for x in st]); tt = 2 * np.array(ls_)
            ax[q].errorbar(tt, m, e, color=col[s], ls=ls, marker=mk, ms=3.5, lw=1.6 if kd == "esc" else 1, capsize=2, label=(f"aspect {s:g}" if s != 1 else "circle") if (kd == "esc" and q == 0) else None)
            L.append(f"aspect {s:4.2f} {kd} {nm:7s}: " + "  ".join(f"{2*l:+d}:{x[0]:+.2f}±{x[1]:.2f}(n{x[2]})" for l, x in zip(ls_, st)))
for q, yl in enumerate(("cos τ (1 = head into the bulk)", "elastic force F₀ (kBT/cell; < 0 restoring)", "cost to move 6 cells in, E₆ (kBT)")):
    ax[q].set_xlabel("time relative to departure (t.u.)"); ax[q].set_ylabel(yl); ax[q].axvline(0, color="0.6", lw=0.8, ls=":")
    for sp_ in ("top", "right"): ax[q].spines[sp_].set_visible(False)
ax[0].legend(fontsize=8, frameon=False, title="solid: escaper, dashed: resident", title_fontsize=8)
n0 = {s: int(((kind == "esc") & (asp == s) & (lag == 0)).sum()) for s in ASP}
fig.suptitle("Depleting escapes in the shape sweep (ellipses + circles, R_eq 1.0/1.5/2.0, 3 seeds): " + ", ".join(f"{s:g}: {n0[s]}" for s in ASP) + " events", fontsize=9.5)
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_etsg_shapes.png"), dpi=200)
# size split of the key contrast at lag -20 frames (-40 t.u.) and 0
L.append("\nsize split (escaper minus resident, same lag):")
for s in ASP:
    for r in sorted(set(req)):
        row = []
        for l in (-20, -10, 0):
            sel = (asp == s) & (req == r) & (lag == l)
            row.append(f"{2*l:+d}: dcos {ms(ct[sel & (kind=='esc')])[0] - ms(ct[sel & (kind=='ctl')])[0]:+.2f} dF0 {ms(F0[sel & (kind=='esc')])[0] - ms(F0[sel & (kind=='ctl')])[0]:+.1f} dE6 {ms(E6[sel & (kind=='esc')])[0] - ms(E6[sel & (kind=='ctl')])[0]:+.0f}")
        L.append(f"aspect {s:4.2f} R_eq {r}: n {int((sel & (kind=='esc')).sum())}  " + " | ".join(row))
open(os.path.join(H, "etsg_shapes.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
