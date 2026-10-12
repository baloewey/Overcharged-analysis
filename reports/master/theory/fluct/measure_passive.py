#!/usr/bin/env python3
"""Passive wall-lattice fluctuations (stage 1 of the elastic + hydrodynamic perturbation theory, 2026-10-06).

For each passive run: +1/2 within BAND cells of the wall, frames with exactly N = 2k of them (all at the wall).
Arclength s (nearest wall point) and depth d. Sites = time-averaged positions, found by sorting each frame in s and choosing the
cyclic label shift that best matches the previous frame (tracking without hops between sites).
Output per run: rms tangential and normal displacement per site (about the site mean), the rms of the collective (mean) tangential
displacement (= rotation of the whole lattice, a zero mode of the continuum theory in a circle), and, for circles, the spectrum
<|u_m|^2> of the tangential displacements, u_m = N^{-1/2} sum_j du_j e^{-2 pi i m j/N}.
Writes fluct_measured.json.
"""
import json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(ROOT, "cavity"))
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline, project  # noqa

BAND, TFRAC = 4.0, 0.2          # band widened for k <= 5 (deeper layer)
RUNS = {**{f"circle_R31.83_k{k:02d}": f"cavity/Rscan_R31.83_passive/k{k:02d}" for k in range(2, 11)},
        "circle_R27.28_k06": "cavity/server_main/main_passive/rho0.07/k06/s01",
        **{f"ellipse_eps{e}_k06": f"cassini/P5_passive_refs/ellipse_eps{e}_k06" for e in ("1.2", "1.5", "2.0")}}


def analyse(run, band=BAND):
    j = json.load(open(os.path.join(run, "input.json"))); bc = j["BC"][0]
    k = bc["customAnch"]; N = 2 * k
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), TFRAC)
    W, s, P = wall_polyline(bc)
    S, D, XY = [], [], []
    for f in fr:
        p = f[f[:, 2] > 0]
        d, sa = project(W, s, p[:, :2])
        m = d < band
        if m.sum() != N or (len(p) != N):
            continue
        o = np.argsort(sa[m]); S.append(sa[m][o]); D.append(d[m][o]); XY.append(p[m][o, :2])
    S, D = np.array(S), np.array(D)
    # tracking: cyclic relabelling to follow the previous frame
    for n in range(1, len(S)):
        best = None
        for sh in range(N):
            cand = np.roll(S[n], -sh)
            ds = (cand - S[n - 1] + P / 2) % P - P / 2
            c = (ds ** 2).sum()
            if best is None or c < best[0]:
                best = (c, sh)
        S[n] = np.roll(S[n], -best[1]); D[n] = np.roll(D[n], -best[1])
    # unwrap per site and remove the site mean
    U = np.unwrap(S * 2 * np.pi / P, axis=0) * P / (2 * np.pi)
    site_s = U.mean(0) % P
    dU = U - U.mean(0)
    dD = D - D.mean(0)
    coll = dU.mean(1)                                           # collective sliding
    dUi = dU - coll[:, None]                                    # internal (relative) displacements
    out = dict(N=N, P=P, frames=len(S), site_s=site_s.tolist(), depth=D.mean(0).tolist(),
               sig_t=dU.std(0).tolist(), sig_t_internal=dUi.std(0).tolist(), sig_n=dD.std(0).tolist(),
               sig_coll=float(coll.std()), drift_total=float(coll[-1] - coll[0]))
    if bc["aInv"][0] == bc["aInv"][1]:
        um = np.fft.fft(dU, axis=1) / np.sqrt(N); vm = np.fft.fft(dD, axis=1) / np.sqrt(N)
        out["spec_t"] = (np.abs(um) ** 2).mean(0)[: N // 2 + 1].tolist()
        out["spec_n"] = (np.abs(vm) ** 2).mean(0)[: N // 2 + 1].tolist()
    # time correlation of the internal tangential modes (for the mobility, stage 3)
    dt = 2.0 * (t[1] - t[0]) / 2.0 if len(t) > 1 else 2.0
    out["dt"] = float(t[1] - t[0])
    ac = []
    for lag in range(0, 60):
        ac.append(float((dUi[: len(dUi) - lag] * dUi[lag:]).mean() / (dUi ** 2).mean()))
    out["acf_internal_t"] = ac
    return out


res = {}
for lab, r in RUNS.items():
    try:
        res[lab] = analyse(os.path.join(ROOT, r), band=9.0 if lab.endswith(('k02', 'k03')) else 6.0 if lab.endswith(('k04', 'k05')) else BAND)
        o = res[lab]
        print(f"{lab:22s} N={o['N']:2d} frames={o['frames']:4d} depth {np.mean(o['depth']):.2f}  sig_t {np.mean(o['sig_t']):.2f} "
              f"(internal {np.mean(o['sig_t_internal']):.2f})  sig_n {np.mean(o['sig_n']):.2f}  collective {o['sig_coll']:.2f}  drift {o['drift_total']:.1f}")
    except Exception as e:
        print(lab, "FAILED", e)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fluct_measured.json"), "w"))
