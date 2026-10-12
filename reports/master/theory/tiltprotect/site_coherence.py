#!/usr/bin/env python3
"""Test B (2026-10-11): do local tilt coherence and registry pinning protect a wall site? Wall split into 2k-2 segments (one per nominal site, arclength), windows of 50 t.u.
(25 frames of 2 t.u., second half). Per segment and window: occupied (a wall +1/2 in the segment in > 80% of frames); for occupied ones: local coherence = mean cos(tau_i - tau_neighbours)
over the defects in the two neighbouring segments (tilt alignment with neighbours), and jitter = rms tangential displacement of the segment's defect along the wall within the window.
Outcome: the segment is EMPTY in the next window (occupied in < 20% of its frames) = local depletion. P(depletion) binned by coherence and by jitter, per class (circle-like /
elongated). Usage: site_coherence.py [run dirs...] (default: the shape sweep) -> site_coherence_<tag>.txt"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; N0 = 2 * k - 2
        if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); spc = P / N0; nf = 25; rows = []
        occ = []; tau = []; spos = []
        for f in fr:
            o = np.zeros(N0, bool); tt = np.full(N0, np.nan); sp = np.full(N0, np.nan)
            if len(f):
                d, j = tree.query(f[:, :2]); w = np.flatnonzero((f[:, 2] > 0) & (d < 4))
                for i in w:
                    seg = int((s[j[i]] + spc / 2) // spc) % N0; o[seg] = True
                    tt[seg] = np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j[i]])))); sp[seg] = s[j[i]]
            occ.append(o); tau.append(tt); spos.append(sp)
        occ, tau, spos = np.array(occ), np.array(tau), np.array(spos)
        for w0 in range(0, len(fr) - 2 * nf, nf):
            A = slice(w0, w0 + nf); B = slice(w0 + nf, w0 + 2 * nf); oa = occ[A].mean(0); ob = occ[B].mean(0)
            for g in range(N0):
                if oa[g] < 0.8: continue
                nb = [(g - 1) % N0, (g + 1) % N0]; c = []
                for m in range(w0, w0 + nf):
                    for h in nb:
                        if np.isfinite(tau[m, g]) and np.isfinite(tau[m, h]): c.append(np.cos(tau[m, g] - tau[m, h]))
                sp = spos[A, g]; sp = sp[np.isfinite(sp)]; sp = (sp - sp[0] + P / 2) % P - P / 2 if len(sp) else sp
                jit = np.std(sp) if len(sp) > 3 else np.nan
                rows.append((np.mean(c) if c else np.nan, jit, float(ob[g] < 0.2)))
        return run, np.array(rows)
    except Exception:
        return run, None
if __name__ == "__main__":
    args = sys.argv[1:]
    if args: runs = [os.path.abspath(a) for a in args]; tag = os.path.basename(os.path.dirname(runs[0]))
    else: runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE")); tag = "shapes"
    man = {os.path.join(SW, r["dir"]): r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=2)
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    groups = {}
    for run, a in res:
        if a is None or not len(a): continue
        if run in man: key = "elongated" if float(man[run]["aspect"]) >= 1.25 else "circle-like"
        else: key = os.path.basename(run)
        groups.setdefault(key, []).append(a)
    for key, arrs in sorted(groups.items()):
        A = np.concatenate(arrs); c, jt, dep = A[:, 0], A[:, 1], A[:, 2]
        say(f"{key}: {len(A)} segment-windows, P(depletion next window) {dep.mean():.4f}")
        say("   by local tilt coherence: " + "  ".join(f"[{a:+.1f},{b:+.1f}): {dep[(c >= a) & (c < b)].mean():.4f} (n {((c >= a) & (c < b)).sum()})" for a, b in ((-1, 0), (0, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 1.01)) if ((c >= a) & (c < b)).sum() > 30))
        say("   by jitter (cells):        " + "  ".join(f"[{a},{b}): {dep[(jt >= a) & (jt < b)].mean():.4f} (n {((jt >= a) & (jt < b)).sum()})" for a, b in ((0, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 99)) if ((jt >= a) & (jt < b)).sum() > 30))
    # joint: coherence x jitter, both classes, and the class ratio within each cell
    if "circle-like" in groups and "elongated" in groups:
        Ac, Ae = np.concatenate(groups["circle-like"]), np.concatenate(groups["elongated"]); CB = (-1, 0.5, 0.8, 1.01); JB = (0, 0.5, 1, 2, 99)
        say("joint (rows: coherence bins -1/0.5/0.8/1; columns: jitter 0/0.5/1/2/inf): P(depletion) circle-like | elongated | ratio")
        tot_w = 0; tot_r = 0
        for a, b in zip(CB[:-1], CB[1:]):
            line = f"   coh [{a:+.1f},{b:+.1f}): "
            for c0, c1 in zip(JB[:-1], JB[1:]):
                mc = (Ac[:, 0] >= a) & (Ac[:, 0] < b) & (Ac[:, 1] >= c0) & (Ac[:, 1] < c1); me = (Ae[:, 0] >= a) & (Ae[:, 0] < b) & (Ae[:, 1] >= c0) & (Ae[:, 1] < c1)
                if mc.sum() > 200 and me.sum() > 200:
                    pc, pe = Ac[mc, 2].mean(), Ae[me, 2].mean(); line += f" {pc:.4f}|{pe:.4f}|{pe / max(pc, 1e-5):4.1f}  "
                    tot_w += me.sum(); tot_r += me.sum() * pc
                else: line += "   (few)              "
            say(line)
        say(f"   elongated P(depletion) observed {Ae[:, 2].mean():.4f}; predicted from circle-like rates at the elongated (coherence, jitter) mix (cells with data) {tot_r / max(tot_w, 1):.4f}; circle-like {Ac[:, 2].mean():.4f}")
    open(os.path.join(H, f"site_coherence_{tag}.txt"), "w").write("\n".join(L) + "\n")
