#!/usr/bin/env python3
"""
Survival of wall (layer) defects and docking of bulk defects, circles vs ovals vs ellipses (master document).

Tracks: +1/2 gap-closed tracks (defect_tracks.build_tracks, rlink 8 cells, gap 4 frames), frames after the first 20%.
State per frame from the distance s to the wall (nearest point of the wall polyline), with hysteresis:
  wall if s < 3 cells, bulk if s > 6 cells, otherwise the previous state (a track starting in between has no state
  until it first crosses one of the thresholds).
Wall visits = maximal runs of wall states. Only visits whose start is observed (the track was in the bulk before) enter
the survival estimate; a visit ends by
  escape   (next state is bulk)                      -> event
  track end before the end of the run (annihilation or tracker loss)  -> censored
  end of the run                                     -> censored
Survival S(tau) of wall visits: Kaplan-Meier estimate with these right-censored durations (pooled over seeds).
Escape rate = escapes per wall defect per time unit (wall defect = time-averaged number of +1/2 in the wall state).
Docking: bulk episodes = maximal runs of bulk states. A bulk episode ends by docking (next state wall), by track end
  (annihilation/loss) or by the end of the run. Reported: docking rate (dockings per time unit per 100 cells of wall),
  docking probability = dockings / (dockings + track ends) over bulk episodes that end inside the run, and the median
  time in the bulk before docking (episodes that started inside the track, i.e. after an escape or creation was observed).
Cached per run in ../data/surv_<name>.npz. Figure ../figs/survival_docking.png.
"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "runs"))
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from defect_tracks import build_tracks          # noqa: E402
from edge_criteria import read_with_angle       # noqa: E402
from wall_fourier import wall_polyline, project  # noqa: E402

DATA = os.path.join(HERE, "..", "data"); FIGS = os.path.join(HERE, "..", "figs")
LA = 22.21; SW, SB = 3.0, 6.0
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
C, K = "cavity/server_main/main_act0.07", "cassini"
GROUPS = {   # family -> list of (label, aspect, [run dirs])
    "circles": [(f"R = {r} ℓa", 1.0, [f"{C}/rho{rho}/k06/s0{s}" for s in (1, 2, 3)])
                for r, rho in (("0.57", "0.15"), ("0.86", "0.10"), ("1.43", "0.06"), ("2.15", "0.04"), ("2.87", "0.03"))],
    "ovals 1.53 ℓa": [(f"b/a {b}", asp, [f"{K}/P3/H34_ba{b}_k06/s0{s}" for s in (1, 2, 3)])
                      for b, asp in (("2.5", 1.18), ("2.0", 1.29), ("1.6", 1.51), ("1.3", 1.94))],
    "ovals 2.16 ℓa": [(f"b/a {b}", asp, [f"{K}/W1_size/H48_ba{b}_k06/s0{s}" for s in (1, 2, 3)])
                      for b, asp in (("2.5", 1.18), ("2.0", 1.29), ("1.6", 1.51), ("1.3", 1.94))],
    "ellipses": [(f"ε {e}", float(e), [f"{K}/W3_ellipse/eps{e}_k06/s0{s}" for s in (1, 2)]) for e in ("1.2", "1.5", "2.0")],
    "lobe oval (k = 3)": [("lobe", 1.38, [f"{K}/P4_lobe_ref/H30.4_ba1.8_k03/s0{s}" for s in (1, 2, 3)])],
}


def analyse(path):
    out = os.path.join(DATA, "surv_" + path.replace("/", "__") + ".npz")
    if os.path.exists(out): return out
    run = os.path.join(RUNS, path)
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); dt = t[1] - t[0]; nfr = len(fr); T = nfr * dt
    W, s, P = wall_polyline(bc)
    tracks, _ = build_tracks(fr, dt, 8.0, 4)
    dur, ev, bulk_dock_t, wall_occ = [], [], [], 0
    nesc = ndock = nb_end = 0
    for tr in tracks.values():
        if tr["q"] < 0: continue
        d, _ = project(W, s, tr["xy"]); n = tr["n"]
        st = np.full(len(d), -1); cur = -1
        for i, di in enumerate(d):
            if di < SW: cur = 1
            elif di > SB: cur = 0
            st[i] = cur
        wall_occ += np.sum(st == 1)
        ends_early = n[-1] < nfr - 1
        i = 0
        while i < len(st):
            if st[i] < 0: i += 1; continue
            j = i
            while j + 1 < len(st) and st[j + 1] == st[i]: j += 1
            observed_start = i > 0 and st[i - 1] >= 0          # entered from the other state
            nxt = st[j + 1] if j + 1 < len(st) else None
            length = (n[j] - n[i] + 1) * dt
            if st[i] == 1:                                      # wall visit
                if nxt == 0: nesc += 1
                if observed_start:
                    dur.append(length); ev.append(1 if nxt == 0 else 0)
            else:                                               # bulk episode
                if nxt == 1:
                    ndock += 1
                    if observed_start or n[i] > 0: bulk_dock_t.append(length)
                elif nxt is None and ends_early:
                    nb_end += 1
            i = j + 1
    nwall = wall_occ / nfr
    np.savez_compressed(out, dur=np.array(dur), ev=np.array(ev), bulk_dock_t=np.array(bulk_dock_t), nesc=nesc, ndock=ndock,
                        nb_end=nb_end, nwall=nwall, T=T, P=P, k=bc.get("customAnch"))
    return out


def km(d, e):
    o = np.argsort(d); d, e = d[o], e[o]
    ts, S, s = [0.0], [1.0], 1.0; nrisk = len(d)
    for tt in np.unique(d):
        m = d == tt; ne = e[m].sum()
        if ne: s *= 1 - ne / nrisk; ts.append(tt); S.append(s)
        nrisk -= m.sum()
    return np.array(ts), np.array(S)


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    paths = [p for g in GROUPS.values() for _, _, ps in g for p in ps]
    with Pool(a.workers) as pool:
        for o in pool.imap_unordered(analyse, paths): print("ok", os.path.basename(o), flush=True)
    load = lambda p: np.load(os.path.join(DATA, "surv_" + p.replace("/", "__") + ".npz"))
    rows = []
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8))
    famcol = {"circles": "#555a62", "ovals 1.53 ℓa": "#2a78d6", "ovals 2.16 ℓa": "#eb6834", "ellipses": "#eda100", "lobe oval (k = 3)": "#1baf7a"}
    fammk = {"circles": "X", "ovals 1.53 ℓa": "o", "ovals 2.16 ℓa": "s", "ellipses": "^", "lobe oval (k = 3)": "D"}
    for fam, items in GROUPS.items():
        for lab, asp, ps in items:
            D = [load(p) for p in ps]
            dur = np.concatenate([d["dur"] for d in D]); ev = np.concatenate([d["ev"] for d in D])
            esc = [float(d["nesc"]) / max(float(d["nwall"]), 1e-9) / float(d["T"]) for d in D]
            dock = [float(d["ndock"]) / float(d["T"]) / float(d["P"]) * 100 for d in D]
            pd = [float(d["ndock"]) / max(float(d["ndock"]) + float(d["nb_end"]), 1) for d in D]
            bt = np.concatenate([d["bulk_dock_t"] for d in D])
            ts, S = km(dur, ev)
            idx = np.searchsorted(ts, 1000, side="right") - 1
            S1000 = S[idx] if len(ts) else np.nan
            half = ts[np.argmax(S <= 0.5)] if (S <= 0.5).any() else np.inf
            rows.append((fam, lab, asp, len(ps), len(dur), int(ev.sum()), np.mean([float(d["nwall"]) for d in D]), np.mean(esc), np.std(esc) / np.sqrt(len(esc)),
                         S1000, half, np.mean(dock), np.std(dock) / np.sqrt(len(dock)), np.mean(pd), np.median(bt) if len(bt) else np.nan))
            panel = {"circles": 0, "ovals 1.53 ℓa": 1, "ellipses": 2, "lobe oval (k = 3)": 2}.get(fam)
            if panel is not None:
                ax = axs[0, panel]
                ax.step(ts, S, where="post", lw=1.6, label=f"{lab} ({len(dur)} visits)" if fam != "lobe oval (k = 3)" else f"lobe-shaped oval, k = 3 ({len(dur)})",
                        color=plt.cm.viridis(len(ax.lines) / 5) if fam != "lobe oval (k = 3)" else "#1baf7a")
    titles = ("(a) circles, k = 6: survival of wall visits", "(b) ovals, half-length 1.53 ℓa", "(c) ellipses (area 2500) and lobe-shaped oval")
    for ax, ti in zip(axs[0], titles):
        ax.set_xscale("log"); ax.set_ylim(0, 1.02); ax.set_xlim(2, 2e4); style(ax)
        ax.set_xlabel("time since the defect reached the wall (time units)", fontsize=9, color=INK)
        ax.set_ylabel("fraction still at the wall, S(τ)", fontsize=9, color=INK); ax.set_title(ti, loc="left", fontsize=10, color=INK)
        ax.legend(frameon=False, fontsize=7)
    for j, (key, ekey, yl, ti) in enumerate(((7, 8, "escapes per wall defect per time unit", "(d) escape rate from the wall"),
                                               (11, 12, "dockings per time unit per 100 cells of wall", "(e) docking rate"),
                                               (13, None, "fraction of bulk +½ episodes ending at the wall", "(f) docking probability (vs annihilation)"))):
        ax = axs[1, j]
        rc = [x for x in rows if x[0] == "circles"]
        for i, x in enumerate(rc):
            xc = 0.93 + 0.035 * i; col = plt.cm.Greys(0.35 + 0.13 * i)
            ax.errorbar([xc], [x[key]], yerr=[x[ekey]] if ekey else None, color=col, marker="X", ms=7, mec="white", lw=0, elinewidth=1.2,
                        capsize=2, label="circles, R = 0.57 → 2.87 ℓa (light → dark)" if i == 0 else None)
            ax.annotate(x[1].replace("R = ", "").replace(" ℓa", ""), (xc, x[key]), xytext=(4, 0), textcoords="offset points", fontsize=6, color=INK2, va="center")
        for fam in [f for f in GROUPS if f != "circles"]:
            r = [x for x in rows if x[0] == fam]
            ax.errorbar([x[2] for x in r], [x[key] for x in r], yerr=[x[ekey] for x in r] if ekey else None, color=famcol[fam],
                        marker=fammk[fam], ms=6, mec="white", lw=1.5, capsize=2, label=fam)
        if key == 7: ax.set_yscale("log")
        ax.set_xlabel("aspect ratio", fontsize=9, color=INK); ax.set_ylabel(yl, fontsize=9, color=INK)
        ax.set_title(ti, loc="left", fontsize=10, color=INK); style(ax)
    axs[1, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Wall-defect survival and docking (+½ tracks, wall: s < 3 cells, bulk: s > 6 cells, hysteresis in between; "
                 "Kaplan–Meier with censored visits; seeds pooled for S, mean ± s.e. over seeds for rates)", fontsize=8.5, color=INK2)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "survival_docking.png"), dpi=170)
    with open(os.path.join(DATA, "survival_docking.txt"), "w") as f:
        f.write("family\tcavity\taspect\tseeds\tvisits\tescapes\twall+1/2\tesc_rate\tesc_se\tS(1000)\tt_half\tdock_rate_per100\tdock_se\tP_dock\tmedian_bulk_t\n")
        for r in rows: f.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    print(open(os.path.join(DATA, "survival_docking.txt")).read())


if __name__ == "__main__":
    main()
