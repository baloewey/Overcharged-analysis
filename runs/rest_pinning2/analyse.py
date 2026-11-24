#!/usr/bin/env python3
"""Passive rest-pinning at crowded spacings (gen.py; user 2026-11-23). Per run, frames from defects.dat (every 2 t.u.), only while the wall still holds the
seeded 2k-2 ring (wall +1/2 count == 2k-2 and no bulk defect yet merged; t < t_break = first time the wall count leaves 2k-2 for >= 10 consecutive frames).
  (1) Tilt relaxation: signed coherent layer tilt arg<exp(i(psi_j - phi_j))> (0 = tails on the homeotropic sites; tracker angle convention of
      rest_pinning/analyse.py) in windows [0,20) [20,60) [60,120) [120,250) [250,500); mean over seeds (circular), with valid-frame counts.
      Question: do seeded 30/60 deg rings relax towards 0 at small spacing (absolute-length pinning) but not at large spacing?
  (2) Per-defect pinning stiffness by equipartition: kappa = kBT / <(tau_j - tau_mean)^2> over valid frames in [20, t_break) of all three tau0
      (per-defect deviation from the ring mean, i.e. the single-defect mode of the mean-field calculation, kBT = 1).
-> summary.txt"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); WIN = [(0, 20), (20, 60), (60, 120), (120, 250), (250, 500)]


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; N = 2 * k - 2
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
    nw = np.array([((f[:, 2] > 0) & (R - np.hypot(f[:, 0], f[:, 1]) < 4)).sum() if len(f) else 0 for f in fr]); bad = nw != N
    tb = t[-1] + 1
    for i in range(len(bad) - 10):
        if bad[i:i + 10].all(): tb = t[i]; break
    Z, dev = {}, []
    for i, f in enumerate(fr):
        if t[i] >= tb or bad[i] or not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        tj = np.angle(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))); zm = np.mean(np.exp(1j * tj))
        for lo, hi in WIN:
            if lo <= t[i] < hi: Z.setdefault((lo, hi), []).append(zm)
        if t[i] >= 20: dev.extend(np.angle(np.exp(1j * tj) / (zm / abs(zm))))
    return dict(k=k, N=N, a=2 * np.pi * R / N, tb=tb, Z=Z, dev=np.array(dev))


if __name__ == "__main__":
    runs = sorted(d for d in glob.glob(H + "/k*_t*_s*") if os.path.exists(d + "/DONE"))
    res = {os.path.basename(r): one(r) for r in runs}
    L = [f"PASSIVE REST PINNING, crowded spacings (R 31.83, no flow; {len(res)} runs finished). Valid frames: wall holds exactly 2k-2, t < t_break.",
         "per (k, tau0): signed coherent tilt by window " + " ".join(f"[{lo},{hi})" for lo, hi in WIN) + " (deg; valid frames) | mean t_break (t.u.)"]
    for k in (4, 6, 9, 10, 12, 14, 16):
        for t0 in (0, 30, 60):
            g = [v for n, v in res.items() if n.startswith(f"k{k:02d}_t{t0:03d}_")]
            if not g: continue
            cells = []
            for wv in WIN:
                zs = [z for v in g for z in v["Z"].get(wv, [])]
                cells.append(f"{np.degrees(np.angle(np.mean(zs))):+6.1f} ({len(zs):3d})" if zs else "    -       ")
            L.append(f"k {k:2d} (spacing {g[0]['a']:4.1f}) tau0 {t0:2d}: " + " | ".join(cells) + f" | t_break {np.mean([v['tb'] for v in g]):6.0f} ({len(g)} seeds)")
    L.append("\nPER-DEFECT PINNING STIFFNESS by equipartition (deviation of each wall defect's tilt from the ring mean, t >= 20 and valid; kBT/rad^2):")
    for k in (4, 6, 9, 10, 12, 14, 16):
        d = np.concatenate([v["dev"] for n, v in res.items() if n.startswith(f"k{k:02d}_")] or [np.array([])])
        if len(d) > 50: L.append(f"  k {k:2d} (spacing {2 * np.pi * 31.830989 / (2 * k - 2):4.1f}): rms deviation {np.degrees(np.sqrt(np.mean(d ** 2))):5.1f} deg -> kappa {1 / np.mean(d ** 2):5.1f} kBT/rad^2 (n {len(d)})")
    open(os.path.join(H, "summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
