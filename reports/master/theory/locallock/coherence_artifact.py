#!/usr/bin/env python3
"""Evidence for the incoherent-frame artefact (user, 2026-10-22). Per director frame (second half) of the active main-sweep runs with director output:
layer coherence C = |<e^{i tau_j}>| over the wall +1/2 (tracker), number of wall +1/2, circulation sign (flow), and the near-wall field
alpha_1 = arg<S e^{2i(theta - phi - pi/2) sgn}>/2 over the cells at u in [d + 0.5, d + 1.5) (d = the run's median wall-defect depth), signed by the
CIRCULATION (meaningful in every frame). Per run: field from all frames (angle of the summed complex average, the local_lock.py estimator) and from
coherent frames only (C > 0.9, exactly 2k-2 wall defects). -> coherence_artifact.npz, coherence_artifact.txt"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; sys.path.insert(0, ROOT + "/runs/cavity"); from edge_criteria import read_with_angle  # noqa
def load(fn, nc):
    with (gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)) as fh: return np.loadtxt((l for l in fh if len(l.split()) == nc and l.split()[0][0].isdigit()))
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        if N < 4: return None
        D = load(run + "/directorfield.dat.gz", 8); F = load(run + "/flowfield.dat.gz", 7)
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        r = np.hypot(X, Y); ph = np.arctan2(Y, X); u = R - r
        deps = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt))); f = fr[i]; rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
            if w.any(): deps.append(np.median(R - rr[w]))
        d = np.median(deps); band = (u >= d + 0.5) & (u < d + 1.5) & (r < R); UB = np.arange(0.0, 6.01, 0.25); rows = []; prof = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            ff = F[np.abs(F[:, 0] - tt) < 1e-6]
            if not len(ff): continue
            xf = ff[:, 1] + 0.5 - Q[0]; yf = ff[:, 2] + 0.5 - Q[1]; rf = np.hypot(xf, yf)
            circ = np.sign(((xf * ff[:, 5] - yf * ff[:, 4]) / np.maximum(rf, 1e-9))[rf < R - 1].mean())
            if circ == 0: continue
            f = fr[i]; rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
            C = abs(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))) if w.any() else 0.0
            dd = D[D[:, 0] == tt]; th = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th[ix, iy] = np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]; sel = band & np.isfinite(th)
            z = np.sum(S[sel] * np.exp(2j * circ * (th[sel] - ph[sel] - np.pi / 2))) / max(sel.sum(), 1)
            rows.append((C, w.sum(), z.real, z.imag))
            pz = []
            for lo_ in UB[:-1]:
                sb = (u >= lo_) & (u < lo_ + 0.25) & (r < R) & np.isfinite(th); pz.append(np.sum(S[sb] * np.exp(2j * circ * (th[sb] - ph[sb] - np.pi / 2))) / max(sb.sum(), 1))
            prof.append(pz)
        if len(rows) < 5: return None
        A = np.array(rows); z = A[:, 2] + 1j * A[:, 3]; coh = (A[:, 0] > 0.9) & (A[:, 1] == N)
        f_all = np.degrees(np.angle(z.mean())); f_coh = np.degrees(np.angle(z[coh].mean())) if coh.sum() >= 3 else np.nan
        P = np.array(prof); return dict(P_all=P.mean(0), P_coh=P[coh].mean(0) if coh.sum() >= 3 else np.full(P.shape[1], np.nan), d=d, run=os.path.relpath(run, ROOT + "/runs"), k=k, a=2 * np.pi * (R - d) / N, C=A[:, 0], nw=A[:, 1], N=N, fz=np.degrees(np.angle(z)), fcoh=coh.mean(), f_all=f_all, f_coh=f_coh)
    except Exception as ex:
        return dict(run=run, err=repr(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s0[1-3]/directorfield.dat.gz"))
    with Pool(8) as p: res = [x for x in p.map(one, runs, chunksize=1) if x and "err" not in x]
    L = [f"{len(res)} runs. Near-wall field 2 alpha_1 (deg) per spacing bin (2 alpha_1 = angle of the complex average, already doubled): fraction of coherent frames | field, all frames | field, coherent frames | field, incoherent frames (C < 0.9)"]
    sp = np.array([x["a"] for x in res])
    for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
        mm = [x for x in res if lo <= x["a"] < hi]
        if not mm: continue
        fz = np.concatenate([x["fz"] for x in mm]); C = np.concatenate([x["C"] for x in mm])
        L.append(f"[{lo:2d},{hi:2d}): {len(mm)} runs, {len(C)} frames | {np.median([x['fcoh'] for x in mm]):.2f} | {np.median([x['f_all'] for x in mm]):5.1f} | {np.nanmedian([x['f_coh'] for x in mm]):5.1f} | {np.median(fz[C < 0.9]) if (C < 0.9).any() else np.nan:5.1f} | depth {np.median([x['d'] for x in mm]):.2f}")
    UC = np.arange(0.125, 6.0, 0.25)
    np.savez("coherence_profiles.npz", uc=UC, a=np.array([x["a"] for x in res]), d=np.array([x["d"] for x in res]), P_all=np.array([x["P_all"] for x in res]), P_coh=np.array([x["P_coh"] for x in res]))
    np.savez("coherence_artifact.npz", a=sp, fcoh=np.array([x["fcoh"] for x in res]), f_all=np.array([x["f_all"] for x in res]), f_coh=np.array([x["f_coh"] for x in res]),
             C=np.concatenate([x["C"] for x in res]), fz=np.concatenate([x["fz"] for x in res]), afr=np.concatenate([np.full(len(x["C"]), x["a"]) for x in res]))
    open("coherence_artifact.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
