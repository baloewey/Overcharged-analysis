#!/usr/bin/env python3
"""Analysis of drive_barrier.csv. U(x) = one-constant cost of dragging a wall +1/2 from its resident depth d0 to depth x (force_vs_N.npz: E_N on delta grid,
circle P 200; N = rho x 200, linear interpolation between grid N). Per run: least-squares fit ln S(x) = a - U(x)/T_a over x from d0 + 0.5 to 10 with
S(x) * n >= 5 (enough events). Classes: coherent (C > 0.85) / incoherent (C < 0.85) from shape_runs.csv; rho levels. -> drive_barrier.txt, fig_drive_barrier.png"""
import csv, numpy as np, matplotlib.pyplot as plt
import os
MODE = os.environ.get("UMODE", "frozen")                                  # frozen (force_vs_N) or relaxed (drag_relaxed: neighbours slide along the wall)
if MODE == "relaxed":
    Z = np.load("../drag_relaxed.npz"); NS = np.array([8, 12, 18, 24, 28]); DEL = Z["DEL"]; EN = np.array([Z[f"relaxed_N{N}"] for N in NS])
else:
    Z = np.load("../force_vs_N.npz"); NS, DEL = Z["NS"], Z["DEL"]; EN = np.array([Z[f"E_N{N}"] for N in NS])
def U(x, d0, rho):
    N = np.clip(rho * 200, NS[0], NS[-1]); E = np.array([np.interp(N, NS, EN[:, j]) for j in range(len(DEL))]); return np.interp(np.clip(x - d0, 0, DEL[-1]), DEL, E)
R = [r for r in csv.DictReader(open("drive_barrier.csv")) if not r.get("err")]; S_ = {r["run"]: r for r in csv.DictReader(open("../../shapes/shape_runs.csv"))}
XG = np.array(sorted(float(k[1:]) for k in R[0] if k.startswith("S"))); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
rows = []
for r in R:
    s = S_[r["run"]]; rho, C, asp, n, d0 = float(s["rho"]), float(s["coh"]), float(s["aspect"]), int(r["n"]), float(r["d0"])
    S = np.array([float(r[f"S{x:g}"]) for x in XG]); ok = (XG >= d0 + 0.5) & (S * n >= 5) & (S > 0)
    Ta = r2 = np.nan
    if ok.sum() >= 3:
        u = np.array([U(x, d0, rho) for x in XG[ok]]); y = np.log(S[ok]); b = np.polyfit(u, y, 1); Ta = -1 / b[0] if b[0] < 0 else np.nan
        r2 = 1 - np.var(y - np.polyval(b, u)) / max(np.var(y), 1e-12)
    rows.append(dict(rho=rho, lev=int(np.argmin(abs(rho - LV))), C=C, asp=asp, d0=d0, S=S, Ta=Ta, r2=r2, n=n))
L = ["DRIVE x BARRIER (excursion depths of wall residents within 20 t.u.; NEW v2 tracker). Medians over runs.",
     "class | rho | runs | d0 | S(d0+1) | S(d0+3) | S(6) | T_a (kBT; fit ln S = a - U/T_a) | R^2 of the fit | runs fitted"]
for cname, sel in (("incoherent C<0.85 (aspect>=1.1)", lambda q: q["C"] < 0.85 and q["asp"] > 1.05), ("coherent C>0.85 (aspect>=1.1)", lambda q: q["C"] > 0.85 and q["asp"] > 1.05), ("circles", lambda q: q["asp"] < 1.05)):
    for i, r0 in enumerate(LV):
        g = [q for q in rows if q["lev"] == i and sel(q)]
        if not g: continue
        s1 = np.median([np.interp(q["d0"] + 1, XG, q["S"]) for q in g]); s3 = np.median([np.interp(q["d0"] + 3, XG, q["S"]) for q in g]); s6 = np.median([np.interp(6, XG, q["S"]) for q in g])
        ta = [q["Ta"] for q in g if np.isfinite(q["Ta"])]; rr = [q["r2"] for q in g if np.isfinite(q["r2"])]
        L.append(f"{cname:32s} | {r0:.2f} | {len(g):3d} | {np.median([q['d0'] for q in g]):.2f} | {s1:.4f} | {s3:.4f} | {s6:.4f} | {np.median(ta) if ta else np.nan:6.1f} | {np.median(rr) if rr else np.nan:.2f} | {len(ta)}")
open(f"drive_barrier_{MODE}.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
fig, ax = plt.subplots(1, 4, figsize=(19, 4.8)); fig.subplots_adjust(wspace=0.32, left=0.05, right=0.99, bottom=0.15, top=0.88)
RC = ("#fdae6b", "#f98c3a", "#e6550d", "#a63603", "#4d1a02")
for j, (cname, sel) in enumerate((("Incoherent ($C<0.85$), aspect $\\geq1.1$", lambda q: q["C"] < 0.85 and q["asp"] > 1.05), ("Coherent ($C>0.85$), all shapes", lambda q: q["C"] > 0.85))):
    for i, r0 in enumerate(LV):
        g = [q for q in rows if q["lev"] == i and sel(q)]
        if len(g) < 3: continue
        Sm = np.median([q["S"] for q in g], 0); d0 = np.median([q["d0"] for q in g])
        ax[j].semilogy(XG, np.where(Sm > 0, Sm, np.nan), "-o", ms=3, color=RC[i], label=fr"$\rho={r0:.2f}$")
        if j == 0:
            u = np.array([U(x, d0, r0) for x in XG]); m = (XG >= d0 + 0.5) & (Sm > 0)
            ax[2].plot(u[m], -np.log(Sm[m]), "-o", ms=3, color=RC[i], label=fr"$\rho={r0:.2f}$")
    ax[j].set_xlabel("Depth $x$ (cells)"); ax[j].set_ylabel("$S(x)$: reach depth $x$ within 20 t.u."); ax[j].set_title(cname, fontsize=11); ax[j].legend(fontsize=8.5)
ax[2].set_xlabel(rf"Elastic cost $U(x)$ to reach depth $x$ ($k_BT$; {MODE} neighbours)"); ax[2].set_ylabel(r"$-\ln S(x)$"); ax[2].set_title("Incoherent layers: collapse on $U$?", fontsize=11); ax[2].legend(fontsize=8.5)
for cname, sel, c, mk in (("Incoherent", lambda q: q["C"] < 0.85 and q["asp"] > 1.05, "#a63603", "o"), ("Coherent", lambda q: q["C"] > 0.85, "#3b6fb6", "s")):
    y = [np.nanmedian([q["Ta"] for q in rows if q["lev"] == i and sel(q)] or [np.nan]) for i in range(5)]; ax[3].plot(LV, y, "-" + mk, color=c, label=cname)
ax[3].set_xlabel(r"Wall density $\rho$"); ax[3].set_ylabel(r"Fitted active temperature $T_a$ ($k_BT$)"); ax[3].set_title("Drive: does $T_a$ depend on $\\rho$?", fontsize=11); ax[3].legend(fontsize=9)
for a, s in zip(ax, "abcd"): a.text(-0.13, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
for e in ("pdf", "png"): fig.savefig(f"fig_drive_barrier_{MODE}.{e}", dpi=170)
