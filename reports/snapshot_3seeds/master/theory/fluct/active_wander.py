#!/usr/bin/env python3
"""How much do wall +1/2 wander in active runs, compared with passive? (sets the realistic disorder level for stage 2)
Per frame: wall +1/2 = within BAND cells of the wall; arclength positions s_j sorted; gaps g_j = s_{j+1} - s_j (cyclic).
  sigma_gap = sd(g)/sqrt(2): for a lattice with independent Gaussian displacements sigma_t about evenly spaced sites, sd(g) = sqrt(2) sigma_t,
  so sigma_gap is directly the sigma_t of the synthetic stage-2 ensembles (no tracking needed); pooled over frames with the modal count.
  Also: distribution of the wall count (vacancies relative to 2k-2 / 2k), depth sd, and for frames with the modal count the tracked
  internal tangential rms (as for passive) and collective rms.
Writes active_wander.txt."""
import glob, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline, project  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
def analyse(run, band=4.0, nsite=None):
    bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.2); W, s, P = wall_polyline(bc)
    cnt, G, D, S = [], [], [], []
    for f in fr:
        p = f[f[:, 2] > 0]; d, sa = project(W, s, p[:, :2]); m = d < band; n = m.sum(); cnt.append(n)
        if n >= 2:
            ss = np.sort(sa[m]); g = np.diff(np.append(ss, ss[0] + P)); G.append((n, g)); D.append((n, d[m])); S.append((n, ss))
    cnt = np.array(cnt); mode = np.bincount(cnt).argmax()
    gm = np.concatenate([g for n, g in G if n == mode]); dm = np.concatenate([d for n, d in D if n == mode])
    sig_gap = gm.std() / np.sqrt(2)
    # tracked internal / collective rms (frames with the modal count)
    Sm = np.array([x for n, x in S if n == mode])
    for i in range(1, len(Sm)):
        best = min(range(mode), key=lambda sh: (((np.roll(Sm[i], -sh) - Sm[i - 1] + P / 2) % P - P / 2) ** 2).sum()); Sm[i] = np.roll(Sm[i], -best)
    U = np.unwrap(Sm * 2 * np.pi / P, axis=0) * P / (2 * np.pi); dU = U - U.mean(0); coll = dU.mean(1)
    return dict(k=k, P=P, mode=mode, frac_mode=np.mean(cnt == mode), cnt_mean=cnt.mean(), cnt_sd=cnt.std(), sig_gap=sig_gap,
                sig_depth=dm.std(), depth=dm.mean(), sig_int=(dU - coll[:, None]).std(), sig_coll=coll.std(), spacing=P / mode)
def line(lab, o):
    return (f"{lab:42s} k={o['k']:2d} wall count {o['cnt_mean']:5.2f}±{o['cnt_sd']:.2f} (mode {o['mode']:2d}, {o['frac_mode']:.0%} of frames), "
            f"spacing {o['spacing']:5.1f}; sigma_gap {o['sig_gap']:.2f}; tracked internal {o['sig_int']:.2f}, collective {o['sig_coll']:.2f}; "
            f"depth {o['depth']:.2f}±{o['sig_depth']:.2f}")
say("== passive references (calibrates sigma_gap against the tracked rms)")
for lab, rel in (("passive circle R31.83 k06", "cavity/Rscan_R31.83_passive/k06"), ("passive circle R27.28 k06", "cavity/server_main/main_passive/rho0.07/k06/s01"),
                 ("passive ellipse eps1.2", "cassini/P5_passive_refs/ellipse_eps1.2_k06"), ("passive ellipse eps1.5", "cassini/P5_passive_refs/ellipse_eps1.5_k06"),
                 ("passive ellipse eps2.0", "cassini/P5_passive_refs/ellipse_eps2.0_k06")):
    say(line(lab, analyse(os.path.join(ROOT, rel))))
say("== active circles, act 0.07, k = 6 (all densities, seeds)")
for f in sorted(glob.glob(ROOT + "/cavity/server_main/main_act0.07/rho*/k06/s[01][0-9]/DONE")):
    run = os.path.dirname(f); R = json.load(open(run + "/input.json"))["BC"][0]["R"]
    try: say(line(f"active circle {run.split('/')[-3]} {run.split('/')[-1]} R {R:.1f}", analyse(run)))
    except Exception as e: say(f"{run}: failed {e}")
say("== active ellipses W3 (area 2500, k = 6)")
for e in ("1.2", "1.5", "2.0"):
    for sd in ("s01", "s02"):
        run = f"{ROOT}/cassini/W3_ellipse/eps{e}_k06/{sd}"
        say(line(f"active ellipse eps{e} {sd}", analyse(run)))
open(os.path.join(H, "active_wander.txt"), "w").write("\n".join(L) + "\n")
