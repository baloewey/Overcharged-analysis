#!/usr/bin/env python3
"""Figures for the tilt-crowding report (2026-10-16). All data from reports/master/theory/{tiltcrowd,shear_tilt,flatpin} outputs and runs/swirl_crowd.
Run: reports/figstyle/run.sh report_figs.py  -> figs/*.pdf/png"""
import os, re, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); T = "/home/baloewe/Physics/Overcharged/reports/master/theory"; F = os.path.join(H, "figs"); LA = 22.21
RED, BLUE, GREEN, GREY = "#c2452d", "#1f5fa8", "#2e7d3c", "0.5"
def save(fig, name): fig.tight_layout(); fig.savefig(os.path.join(F, name + ".png"), dpi=200); plt.close(fig)
def clean(ax): [ax.spines[s].set_visible(False) for s in ("top", "right")]
def logx(ax, ticks=(5, 7, 10, 15, 20, 30, 50)):
    ax.set_xscale("log"); ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.minorticks_off()
def rows(fn, prefix):
    return [[c.strip() for c in l.split("|")] for l in open(fn) if l.startswith(prefix)]

# ---------- data ----------
A = np.load(f"{T}/tiltcrowd/tilt_crowding_runs.npz"); ak, aR, arho, atc, ata = A["k"], A["R"], A["rho"], A["tilt_coh"], A["tilt_all"]; aa = 2 * np.pi * aR / (2 * ak - 2)
P = np.load(f"{T}/tiltcrowd/passive_runs.npz")
S = np.load(f"{T}/tiltcrowd/shape_tilt.npz", allow_pickle=True)
DFm = rows(f"{T}/tiltcrowd/defect_flow_main10.txt", "rho")
WSm = rows(f"{T}/tiltcrowd/wallshear_main10.txt", "rho")

# ---------- Fig. 1: the observation ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
pts = {}
for k_, R_, t_ in zip(ak, aR, atc):
    if np.isfinite(t_): pts.setdefault((k_, round(R_, 2)), []).append(t_)
for (k_, R_), v in pts.items():
    ax[0].errorbar(2 * np.pi * R_ / (2 * k_ - 2), np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0, fmt="o", ms=3.5, color=plt.cm.plasma(min((R_ / LA - 0.8) / 1.5, 1) * 0.8), capsize=1.5)
ax[0].axhline(60, color="k", ls="--", lw=0.8); ax[0].text(48, 61, r"arccos$(1/\lambda)$", ha="right", fontsize=8)
logx(ax[0]); ax[0].set_xlabel("wall-defect spacing $a=2\\pi R/(2k-2)$ (cells)"); ax[0].set_ylabel("tilt of the coherent wall layer (deg)")
ax[0].set_title("(a) Active circles, 10 seeds (colour: R/$\\ell_a$ 0.8 → 2.2)", fontsize=10)
trip, idx, fam, asp = S["trip"], S["run_idx"], S["fam"], S["asp"]
for lab, m, col, mk in (("circles (shape sweep)", fam == "circle", "k", "o"), ("ellipses, aspect ≥ 1.5", (fam == "ellipse") & (asp >= 1.5), RED, "s"), ("Cassini ovals, aspect ≥ 1.5", (fam == "oval") & (asp >= 1.5), BLUE, "^")):
    xs, ys = [], []
    for i in np.flatnonzero(m):
        t_ = trip[idx == i]
        if len(t_) > 500: xs.append(np.median(t_[:, 0])); ys.append(np.degrees(np.median(t_[:, 2])))
    ax[1].plot(xs, ys, mk, color=col, ms=3.5, alpha=0.7, label=lab)
ax[1].axhline(60, color="k", ls="--", lw=0.8); logx(ax[1]); ax[1].set_xlabel("local site spacing $a(s)$, run median (cells)"); ax[1].set_ylabel("local tilt (deg, run median)")
ax[1].set_title("(b) Elongated cavities (3 seeds)", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
for Aamp, mk, lab in ((0.01, "s", "passive, imposed vortex 0.01"), (0.03, "D", "passive, imposed vortex 0.03")):
    m = (P["A"] == Aamp) & np.isfinite(P["tilt"]); ax[2].plot(P["a"][m], P["tilt"][m], mk, color=BLUE if Aamp == 0.01 else GREEN, ms=5, label=lab)
m31 = np.abs(aR - 31.83) < 0.05
for k_ in sorted(set(ak[m31])):
    pass
act = {}
for r in rows(f"{T}/tiltcrowd/bulk_spiral.txt", " ") + rows(f"{T}/tiltcrowd/bulk_spiral.txt", "1") + rows(f"{T}/tiltcrowd/bulk_spiral.txt", "2"):
    try: act[int(r[0])] = (float(r[2]), float(r[3]))
    except Exception: pass
ka = sorted(act); ax[2].plot([act[k][0] for k in ka], [act[k][1] for k in ka], "o-", color=RED, ms=4, label="active, same R (10 seeds)")
ax[2].axhline(60, color="k", ls="--", lw=0.8); logx(ax[2]); ax[2].set_xlabel("wall-defect spacing (cells)"); ax[2].set_ylabel("tilt of the coherent wall layer (deg)")
ax[2].set_title("(c) Passive layers drop too (R = 1.43 $\\ell_a$)", fontsize=10); ax[2].legend(fontsize=8, frameon=False)
for a_ in ax: clean(a_)
save(fig, "fig1_observation")

# ---------- Fig. 2: measurement ----------
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
FT = [r for r in csv.DictReader(open("/home/baloewe/Physics/Overcharged/runs/cavity/analysis_flow/flow_tilt_3seeds.csv"))]
g = lambda r, c: float(r[c]) if r[c] not in ("", "nan") else np.nan
rb = [(0.025, 0.055), (0.055, 0.075), (0.075, 0.095), (0.095, 0.11), (0.11, 0.13), (0.13, 0.16)]; xm = [0.04, 0.065, 0.085, 0.10, 0.12, 0.145]
win = (ak >= 2) & (ak <= 10) & (aR / LA >= 0.8) & (aR / LA <= 2.2)
ax[0].plot(xm, [np.nanmedian(atc[win & (arho >= lo) & (arho < hi)]) for lo, hi in rb], "o-", color=RED, label="coherent frames only ($C>0.9$), 10 seeds")
ax[0].plot(xm, [np.nanmedian(ata[win & (arho >= lo) & (arho < hi)]) for lo, hi in rb], "s-", color=GREY, label="all frames, |frame mean|, 10 seeds")
for col_, lab, mk, cc in (("tilt_abs", "old table: mean |τ| over defects", "^", BLUE), ("tilt_flow", "old table: frame mean signed by circulation", "v", GREEN)):
    vals = []
    for lo, hi in rb:
        v = [np.degrees(g(r, col_)) for r in FT if 2 <= g(r, "k") <= 10 and 0.8 <= g(r, "R") / LA <= 2.2 and lo <= g(r, "rho") < hi]
        vals.append(np.nanmedian(v) if v else np.nan)
    ax[0].plot(xm, vals, mk + "--", color=cc, label=lab + " (3 seeds)")
ax[0].set_xlabel("wall density $\\rho=k/\\pi R$"); ax[0].set_ylabel("tilt (deg)"); ax[0].legend(fontsize=7.5, frameon=False); ax[0].set_title("(a) Estimators of the tilt", fontsize=10)
LO = rows(f"{T}/tiltcrowd/lock_offset.txt", "[")
mids = [np.sqrt(float(re.findall(r"[\d.]+", r[0])[0]) * float(re.findall(r"[\d.]+", r[0])[1])) for r in LO]
ax[1].plot(mids, [float(r[3]) for r in LO], "o-", color=RED, label="orientation offset from the lock, c")
ax2 = ax[1].twinx(); ax2.plot(mids, [float(r[2]) for r in LO], "s--", color=GREY, label="tilt"); ax2.set_ylabel("tilt (deg)", color=GREY)
logx(ax[1]); ax[1].axhline(0, color="0.8"); ax[1].set_xlabel("wall-defect spacing (cells)"); ax[1].set_ylabel("lock offset c (deg, signed with the tilt)", color=RED)
ax[1].set_title("(b) Cores are nearly relaxed on their positions", fontsize=10)
for a_ in ax: clean(a_)
save(fig, "fig2_measurement")

# ---------- Fig. 3: flow ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
W31 = [r for r in WSm if abs(float(r[2]) - 31.8) < 0.1]
ax[0].plot([float(r[3]) for r in W31], [float(r[4]) for r in W31], "o", color=RED, ms=3.5, label="active: interior $|\\langle v_\\theta\\rangle|$")
ax[0].plot([float(r[3]) for r in W31], [float(r[7]) for r in W31], "s", color=GREY, ms=3.5, label="active: strain/rotation in the wall band")
ax[0].plot(P["a"], P["vth"] / np.nanmax(P["vth"]) * 0.0 + P["vth"], "D", color=BLUE, ms=4, label="passive: interior $|\\langle v_\\theta\\rangle|$")
logx(ax[0]); ax[0].set_xlabel("wall-defect spacing (cells)"); ax[0].set_ylabel("flow measure"); ax[0].set_title("(a) Mean flow does not change (R = 1.43 $\\ell_a$)", fontsize=10); ax[0].legend(fontsize=7.5, frameon=False)
da = np.array([float(r[2]) for r in DFm]); dphi2 = np.array([float(r[5].split()[1]) for r in DFm]); dexx = np.array([float(r[4].split()[2]) for r in DFm]); dtau2 = np.array([float(r[5].split()[2]) for r in DFm]); dT = np.array([float(r[6]) for r in DFm])
ax[1].plot(da, dphi2, "o", color=RED, ms=2.5, alpha=0.5, label="active (10-seed window)"); ax[1].plot(P["a"], P["phi2"], "D", color=BLUE, ms=4.5, label="passive")
logx(ax[1]); ax[1].axhline(0, color="0.8"); ax[1].set_xlabel("wall-defect spacing (cells)"); ax[1].set_ylabel("rotation of the strain axes, 2 cells inside (deg)")
ax[1].set_title("(b) Local flow next to the layer", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
ok = np.isfinite(dT)
ax[2].plot(dtau2[ok], dT[ok], "o", color=RED, ms=2.5, alpha=0.5, label="active"); ax[2].plot(P["tau2"], P["tilt"], "D", color=BLUE, ms=4.5, label="passive")
ax[2].plot([35, 80], [35, 80], "k--", lw=0.8); ax[2].set_xlabel("Leslie angle of the local flow, 2 cells inside (deg)"); ax[2].set_ylabel("measured tilt (deg)")
ax[2].set_title("(c) Local Leslie angle vs measured tilt", fontsize=10); ax[2].legend(fontsize=8, frameon=False)
for a_ in ax: clean(a_)
save(fig, "fig3_flow")

# ---------- Fig. 4: elasticity ----------
fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))
FA = [l.split("|") for l in open(f"{T}/flatpin/flat_pinning_amp.txt") if re.match(r"\s*\d", l)]
aF = [float(r[0]) for r in FA]
for j, (d, c) in enumerate(((1.0, "#9ecae1"), (1.5, "#4292c6"), (2.0, "#08519c"))):
    ax[0].plot(aF, [float(r[j + 1]) for r in FA], "o-", color=c, label=f"depth d = {d}")
ax[0].axhline(0, color="0.7"); ax[0].set_xlabel("spacing a (cells)"); ax[0].set_ylabel("pinning amplitude A (kBT per defect)")
ax[0].set_title("(a) Flat wall, first order in $K_1-K_3$\n$E=A\\cos\\tau$: A > 0 prefers τ = 180°", fontsize=9.5); ax[0].legend(fontsize=8, frameon=False)
CE = [[c.strip() for c in l.split("|")] for l in open(f"{T}/tiltcrowd/crowd_elastic.txt") if re.match(r"\s*\d", l)]
ce = [(float(r[1]), float(r[2].split()[0]), float(r[3].split()[0].split("/")[0])) for r in CE if len(r) > 3]
ax[1].plot([c[0] for c in ce], [c[1] for c in ce], "o-", color=RED, label="registry (ΔK pinning)"); ax[1].plot([c[0] for c in ce], [c[2] for c in ce], "s-", color=BLUE, label="core twist")
logx(ax[1]); ax[1].set_xlabel("wall-defect spacing (cells)"); ax[1].set_ylabel("stiffness per defect (kBT/rad$^2$)"); ax[1].set_title("(b) Circle (R 31.8): both soften with crowding", fontsize=9.5); ax[1].legend(fontsize=8, frameon=False)
FD = [l for l in open(f"{T}/shear_tilt/flatwall_depth.txt") if l.startswith("K1/K3")]
for ratio, mk in (("1.00", "o"), ("1.82", "s")):
    for s_, c in (("0.003", "#9ecae1"), ("0.010", "#4292c6"), ("0.030", "#08519c"), ("0.100", "k")):
        L_ = [l for l in FD if f"K1/K3 {ratio} s {s_}" in l]
        xs = [float(re.search(r"a_s\s+(\d+)", l).group(1)) for l in L_]; ys = [float(re.findall(r"tilt\s+([-\d.]+)", l)[-1]) for l in L_]
        ax[2].plot(xs, ys, mk + ("-" if ratio == "1.00" else "--"), color=c, ms=4, label=f"K1/K3 {ratio}, s {s_}" if s_ in ("0.003", "0.100") else None)
ax[2].axhline(60, color="k", ls=":", lw=0.8); ax[2].set_xlabel("site spacing $a_s$ (cells)"); ax[2].set_ylabel("tilt (deg)"); ax[2].set_title("(c) Flat-wall Q-tensor in shear (converged)\nsolid: one-constant; dashed: K1/K3 = 1.82", fontsize=9.5); ax[2].legend(fontsize=7, frameon=False)
B = np.load(f"{T}/tiltcrowd/bulk_spiral_runs.npz")
sc = ax[3].scatter(B["lean85"], B["tilt"], c=np.log(B["a"]), cmap="viridis", s=12); ax[3].set_xlabel("ring lean at 0.85 R (deg past tangential)"); ax[3].set_ylabel("wall-layer tilt (deg)")
ax[3].set_title("(d) Ring lean is coupled to the tilt\n(colour: log spacing; R 31.8, 116 runs)", fontsize=9.5)
for a_ in ax: clean(a_)
save(fig, "fig4_elastic")

# ---------- Fig. 5: depth ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
dd = np.array([float(r[3]) for r in DFm]); ws = {r[0]: float(r[5]) + float(r[6]) for r in WSm}; ds = np.array([ws.get(r[0], np.nan) for r in DFm])
ok = np.isfinite(dT) & np.isfinite(ds)
ax[0].plot(da, dd, "o", color=RED, ms=2.5, alpha=0.5, label="active"); ax[0].plot(P["a"], P["d"], "D", color=BLUE, ms=4.5, label="passive")
logx(ax[0]); ax[0].set_xlabel("wall-defect spacing (cells)"); ax[0].set_ylabel("core depth d below the wall (cells)"); ax[0].set_title("(a) Crowded cores sit closer to the wall", fontsize=10); ax[0].legend(fontsize=8, frameon=False)
for lab, x, y, c, mk in (("active", dd[ok] * np.sqrt(ds[ok]), dT[ok], RED, "o"), ("passive", P["d"] * np.sqrt(P["s"]), P["tilt"], BLUE, "D")):
    cf = np.polyfit(x, y, 1); ax[1].plot(x, y, mk, color=c, ms=3 if lab == "active" else 4.5, alpha=0.5 if lab == "active" else 1, label=f"{lab}: tilt = {cf[1]:.0f} + {cf[0]:.0f} d√s")
    xx = np.linspace(x.min(), x.max(), 10); ax[1].plot(xx, np.polyval(cf, xx), color=c, lw=1)
    ax[2].plot(x, y - 60, mk, color=c, ms=3 if lab == "active" else 4.5, alpha=0.5 if lab == "active" else 1, label=lab)
ax[1].axhline(60, color="k", ls="--", lw=0.8); ax[1].set_xlabel("$d\\sqrt{s}$ (∝ depth / alignment length)"); ax[1].set_ylabel("tilt (deg)"); ax[1].set_title("(b) One variable orders both sets", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
ax[2].axhline(0, color="k", lw=0.8); ax[2].set_xlabel("$d\\sqrt{s}$"); ax[2].set_ylabel("tilt − arccos(1/λ) (deg)")
ax[2].set_title("(c) The correction term can only be ≥ 0;\ncrowded layers sit below", fontsize=10); ax[2].legend(fontsize=8, frameon=False)
for a_ in ax: clean(a_)
save(fig, "fig5_depth")
print("done")

# ---------- Fig. 6: continuum vs MPCD across spacing ----------
def parse_disk(fn, model=None, U=None):
    out = []
    for l in open(fn):
        m = re.search(r"(\w+)\s+k\s+(\d+) \(spacing\s+([\d.]+)\) U ([\d.]+): .*tilt\s+([-+\d.]+)", l)
        if m and (model is None or m.group(1) == model) and (U is None or float(m.group(4)) == U): out.append((float(m.group(3)), float(m.group(5))))
    return np.array(out)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot([act[k][0] for k in ka], [act[k][1] for k in ka], "o-", color=RED, lw=2, label="MPCD active (10 seeds)")
m = (P["A"] == 0.03) & np.isfinite(P["tilt"]); ax[0].plot(P["a"][m], P["tilt"][m], "D", color=GREEN, ms=5, label="MPCD passive, imposed vortex 0.03")
for fn, mod, lab, c, ls in ((f"{T}/shear_tilt/disk_spacing.txt", "oneK", "continuum disk, one constant", "0.3", "--"),
                            (f"{T}/shear_tilt/disk_spacing.txt", "aniso_mpcd", "continuum disk, $K_1/K_3=1.82$", BLUE, "--"),
                            (f"{T}/shear_tilt/disk_spacing_055.txt", "aniso_mpcd", "continuum disk, $K_1/K_3=0.55$ (pins at τ = 0)", "#8a5ac2", "-")):
    d = parse_disk(fn, mod, 1.0); ax[0].plot(d[:, 0], d[:, 1], "s" + ls, color=c, ms=4, label=lab)
AS = [(float(re.search(r"spacing\s+([\d.]+)", l).group(1)), float(re.search(r"tilt\s+([-+\d.]+)", l).group(1)), float(re.search(r"zeta ([\d.]+)", l).group(1))) for l in open(f"{T}/shear_tilt/active_spacing.txt") if l.startswith("zeta")]
for z, mk in ((0.3, "^"), (1.0, "v")):
    d = np.array([(a_, t_) for a_, t_, zz in AS if zz == z]); ax[0].plot(d[:, 0], d[:, 1], mk + ":", color="k", ms=4, label=f"active continuum, one constant, ζ {z}")
ax[0].axhline(60, color="k", ls="--", lw=0.6); logx(ax[0]); ax[0].set_xlabel("wall-defect spacing (cells)"); ax[0].set_ylabel("tilt of the wall layer (deg)")
ax[0].set_title("(a) Circles, R = 31.8: MPCD vs continuum", fontsize=10); ax[0].legend(fontsize=7, frameon=False, loc="lower right")
for fn, rt, lab, c in ((f"{T}/shear_tilt/flatwall_depth.txt", "1.00", "one constant, full wall order", "0.3"), (f"{T}/shear_tilt/flatwall_depth.txt", "1.82", "$K_1/K_3=1.82$", BLUE),
                       (f"{T}/shear_tilt/flatwall_depth_melt055.txt", "0.55", "$K_1/K_3=0.55$, full order", "#8a5ac2")):
    xs, ys = [], []
    for l in open(fn):
        if l.startswith(f"K1/K3 {rt}") and ("s 0.030" in l) and ("S_w 0.60" not in l and "S_w 0.40" not in l and "S_w 0.20" not in l):
            xs.append(float(re.search(r"a_s\s+(\d+)", l).group(1))); ys.append(float(re.findall(r"tilt\s+([-+\d.na]+)", l)[-1]))
    ax[1].plot(xs, ys, "o-", color=c, label=lab)
xs, ys = [], []
for l in open(f"{T}/shear_tilt/flatwall_depth_melt055.txt"):
    if "S_w 0.60" in l: xs.append(float(re.search(r"a_s\s+(\d+)", l).group(1))); ys.append(float(re.findall(r"tilt\s+([-+\d.na]+)", l)[-1]))
ax[1].plot(xs, ys, "o--", color="#8a5ac2", label="$K_1/K_3=0.55$, wall order 0.6")
ax[1].axhline(60, color="k", ls="--", lw=0.6); ax[1].set_xlabel("site spacing $a_s$ (cells)"); ax[1].set_ylabel("tilt (deg)"); ax[1].set_title("(b) Flat wall in shear: the pinning sign decides", fontsize=10)
ax[1].legend(fontsize=7.5, frameon=False)
for a_ in ax: clean(a_)
save(fig, "fig6_continuum")

# ---------- Fig. 7: tilt potential across crowding + band profile ----------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
KS = (3, 4, 5, 6, 8, 10, 12, 15, 20)
def ks_rows(fn):
    out = {}; key = None
    for l in open(fn):
        if l.startswith("-- bulk"): key = l.strip()[3:]; out[key] = []
        m = re.search(r"k\s+(\d+) \(spacing\s+([\d.]+)\): A1\s+([-+\d.]+).*?A2\s+([-+\d.]+)", l)
        if m and key: out[key].append((int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return out
KR = ks_rows(f"{T}/aniso2/k_scan2.txt")
for key, c in zip(sorted(KR), (RED, "#e08a3c", BLUE, "#6baed6", "k")):
    d = np.array(KR[key]); N = 2 * d[:, 0] - 2; ax[0].plot(d[:, 1], d[:, 2] / N, "o-", color=c, ms=3.5, label=key)
    ax[1].plot(d[:, 1], d[:, 3] / N, "o-", color=c, ms=3.5, label=key)
for a_, yl, tt in ((ax[0], "$A_1$ per wall defect ($k_BT$)", "(a) First harmonic ($A_1<0$: pull towards τ = 0)"), (ax[1], "$A_2$ per wall defect ($k_BT$)", "(b) Second harmonic (needed for a finite tilt): not robust")):
    a_.axhline(0, color="0.7"); logx(a_); a_.set_xlabel("wall-defect spacing (cells)"); a_.set_ylabel(yl); a_.set_title(tt, fontsize=10)
ax[0].legend(fontsize=6.5, frameon=False)
BP = np.load(f"{T}/tiltcrowd/band_profile.npz"); cmap = plt.cm.viridis(np.linspace(0, 0.9, len(BP["ks"])))
for i, (k_, c) in enumerate(zip(BP["ks"], cmap)):
    if k_ in (3, 6, 10, 14, 20): ax[2].plot(BP["mid"], np.where(BP["unif"][i] > 0.3, BP["alpha"][i], np.nan), "o-", color=c, ms=3, label=f"k = {k_}")
ax[2].axhline(30, color="k", ls="--", lw=0.8); ax[2].text(14.5, 30.6, "Leslie angle θ$_L$", ha="right", fontsize=8)
ax[2].set_xlabel("depth below the wall (cells)"); ax[2].set_ylabel("director angle from the tangent (deg)"); ax[2].set_title("(c) Director across the wall band (MPCD, R 31.8)", fontsize=10); ax[2].legend(fontsize=7.5, frameon=False)
for a_ in ax: clean(a_)
save(fig, "fig7_potential_band")
print("done 6-7")
