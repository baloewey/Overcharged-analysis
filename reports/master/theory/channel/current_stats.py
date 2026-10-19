#!/usr/bin/env python3
"""Topological current vs channel width, with statistics (user, 2026-10-15). All active winding channels at L 192 with wall charges +6/-6 (mixed) or +6/+6
(symmetric control): runs/channel_current (6000 t.u.), channel_mixedW, channel_mixed (W 40), channel_narrow and channel (controls). Per run: per-line current
J = net mid-line crossings per unit length per 1000 t.u. (second half), and the motion/event decomposition (current_decomp.one; integrated terms / bulk width).
Per width and type: mean, s.e.m., n; mixed minus control with a t-statistic. -> current_stats.txt, fig_current_stats.png (+pdf)"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); R_ = "/home/baloewe/Physics/Overcharged/runs"
from current_decomp import one  # noqa


def classify(run):
    J = json.load(open(run + "/input.json")); BC = J["BC"]; W = BC[3]["Q"][1] - BC[2]["Q"][1]; act = J["species"][0]["act"]
    cb, ctop = BC[2].get("customAnch", 0), -BC[3].get("customAnch", 0)      # wall charges: bottom +customAnch, top -customAnch
    if act == 0 or J["domain"][0] != 192: return None
    if abs(J["species"][0]["pop"] / (J["domain"][0] * (W + 2)) - 20) > 1: return None     # wrong particle number (channel_current W16/W32, 2026-10-16)
    if cb == 6 and ctop == -6: return ("mixed", W)
    if cb == 6 and ctop == 6: return ("sym", W)
    return None


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cands = sorted(set(glob.glob(f"{R_}/channel_current/*/") + glob.glob(f"{R_}/channel_fix/W*/") + glob.glob(f"{R_}/channel_curve/W*/") + glob.glob(f"{R_}/channel_mixedW/*/") + glob.glob(f"{R_}/channel_mixed/mixed0*/") +
                       glob.glob(f"{R_}/channel_narrow/W*_m6_*/") + glob.glob(f"{R_}/channel/W*_m6_L192_*/")))
    runs = [r.rstrip("/") for r in cands if os.path.exists(r + "DONE")]; cls = {r: classify(r) for r in runs}; runs = [r for r in runs if cls[r] and cls[r][1] >= 14]
    with Pool(int(os.environ.get("NPROC", 8))) as pool: res = dict(zip(runs, pool.map(one, runs)))
    G = {}
    for r in runs:
        typ, W = cls[r]; x = res[r]; wb = W - 8; G.setdefault((typ, W), []).append((x["cross"], (sum(x["Jm"].values()) + x["ev_c"] + x["ev_a"]) / wb,
            (x["Jm"][(1, "wall")]) / wb, (x["Jm"][(1, "pair")] + x["Jm"][(1, "other")]) / wb, (x["Jm"][(-1, "pair")] + x["Jm"][(-1, "other")] + x["Jm"][(-1, "wall")]) / wb, (x["ev_c"] + x["ev_a"]) / wb))
    L = ["TOPOLOGICAL CURRENT vs width (L 192; per-line current per unit length per 1000 t.u.; mean ± s.e.m. over seeds)",
         "W | type | n | crossings | motion+events | of which: wall-ejected +1/2, other +1/2, -1/2, pair events"]
    Ws = sorted(set(W for _, W in G))
    for W in Ws:
        for typ in ("mixed", "sym"):
            if (typ, W) not in G: continue
            A = np.array(G[(typ, W)]); m = A.mean(0); se = A.std(0, ddof=1) / np.sqrt(len(A)) if len(A) > 1 else np.full(A.shape[1], np.nan)
            L.append(f"{W:3.0f} | {typ:5s} | {len(A):2d} | {m[0]:+.4f} ± {se[0]:.4f} | {m[1]:+.4f} ± {se[1]:.4f} | {m[2]:+.4f} {m[3]:+.4f} {m[4]:+.4f} {m[5]:+.4f}")
        if ("mixed", W) in G and ("sym", W) in G:
            a, b = np.array(G[("mixed", W)])[:, 0], np.array(G[("sym", W)])[:, 0]
            if len(a) > 1 and len(b) > 1:
                tstat = (a.mean() - b.mean()) / np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)); L.append(f"      mixed − control: {a.mean()-b.mean():+.4f}, t = {tstat:.1f}")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for typ, col, lab in (("mixed", "#c2452d", "mixed walls (+6/$-$6)"), ("sym", "0.45", "symmetric control (+6/+6)")):
        ws = [W for W in Ws if (typ, W) in G]; A = [np.array(G[(typ, W)])[:, 0] for W in ws]
        ax[0].errorbar(ws, [a.mean() for a in A], [a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0 for a in A], fmt="o-", color=col, capsize=3, label=lab)
        for W, a in zip(ws, A): ax[0].plot([W] * len(a), a, ".", color=col, alpha=0.4)
    ax[0].axhline(0, color="0.8", lw=0.8); ax[0].set_xlabel("channel width W (cells)"); ax[0].set_ylabel("current across the mid-line\n(per unit length per 1000 t.u.)")
    ax[0].set_title("(a) Topological current vs width", fontsize=10); ax[0].legend(fontsize=8, frameon=False)
    ws = [W for W in Ws if ("mixed", W) in G]; M = np.array([np.array(G[("mixed", W)]).mean(0) for W in ws])
    for j, (lab, col) in enumerate((("$+1/2$ ejected from the + wall", "#c2452d"), ("other $+1/2$", "#e08a3c"), ("$-1/2$", "#1f5fa8"), ("pair creation/annihilation", "0.4")), start=2):
        ax[1].plot(ws, M[:, j], "o-", color=col, label=lab)
    ax[1].plot(ws, M[:, 1], "k--", lw=1, label="total"); ax[1].axhline(0, color="0.8", lw=0.8)
    ax[1].set_xlabel("channel width W (cells)"); ax[1].set_ylabel("contribution to the current"); ax[1].set_title("(b) Mixed channels: who carries the current", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_current_stats.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "current_stats.txt"), "w").write("\n".join(L) + "\n")
