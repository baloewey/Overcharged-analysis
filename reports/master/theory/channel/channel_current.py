#!/usr/bin/env python3
"""Topological current and defect cycling in winding channels (user, 2026-10-13: 'topological capacitor'). Works on any eMode-4 channel run (also homeotropic).
Per run (second half, frames every 2 t.u.):
 * wall charges imposed: bottom +c_b, top -c_t (from customAnch); layer occupancy: +1/2 and -1/2 within 4 cells of each wall.
 * charge profile across the width (10 bins): net charge, +1/2, -1/2 per frame.
 * CURRENT of topological charge in +y: (a) mid-line crossings of linked defects (same sign, < 2.5 cells, periodic in x): J_cross = sum q * (up - down) per unit length per 1000 t.u.;
   (b) J_flux = <sum_bulk q_i v_y,i> / L (v from the displacement over one frame), per 1000 t.u.
 * pair NUCLEATION (new +1/2 and new -1/2 within 3 cells) and ANNIHILATION (a +1/2 and a -1/2 within 3 cells both vanish) events: rates and their positions across the width
   (fraction within 4 cells of the bottom wall / bulk / top wall). Also 'layer recombination': a vanishing bulk defect next to a vanishing opposite-sign layer defect.
Usage: channel_current.py RUNDIR [...] -> prints one block per run; with --group, a summary grouped by directory pattern (W, charges, act)."""
import sys, os, json, re, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
def analyse(run):
    J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; W = y1 - y0
    cb, ct = J["BC"][2].get("customAnch", 0), J["BC"][3].get("customAnch", 0); Qbw, Qtw = cb, -ct; act = J["species"][0]["act"]
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); t = np.array(t); dt = float(np.median(np.diff(t))); T = len(fr) * dt; ym = 0.5 * (y0 + y1)
    nb = 10; edges = np.linspace(0, 1, nb + 1); hp = np.zeros(nb); hm = np.zeros(nb); occ = np.zeros(4); nf = 0
    up = down = 0.0; flux = []; nuc_y, ann_y = [], []
    def link(f, g):                       # nearest same-sign partner in g for each defect of f (periodic in x); returns index or -1
        out = -np.ones(len(f), int)
        for i in range(len(f)):
            s = np.flatnonzero(np.sign(g[:, 2]) == np.sign(f[i, 2]))
            if not len(s): continue
            dx = (g[s, 0] - f[i, 0] + L / 2) % L - L / 2; dd = np.hypot(dx, g[s, 1] - f[i, 1]); k = np.argmin(dd)
            if dd[k] < 2.5: out[i] = s[k]
        return out
    prev = None
    for f in fr:
        if not len(f): prev = None; continue
        nf += 1; u = (f[:, 1] - y0) / W; hp += np.histogram(u[f[:, 2] > 0], edges)[0]; hm += np.histogram(u[f[:, 2] < 0], edges)[0]
        nearb = f[:, 1] - y0 < 4; neart = y1 - f[:, 1] < 4
        occ += [((f[:, 2] > 0) & nearb).sum(), ((f[:, 2] < 0) & nearb).sum(), ((f[:, 2] > 0) & neart).sum(), ((f[:, 2] < 0) & neart).sum()]
        if prev is not None and len(prev):
            fw = link(prev, f); bw = link(f, prev); s_ = 0.0
            for i, j in enumerate(fw):
                if j < 0: continue
                q = 0.5 * np.sign(prev[i, 2]); ya, yb = prev[i, 1], f[j, 1]
                if ya < ym <= yb: up += q
                elif yb < ym <= ya: down += q
                if not (prev[i, 1] - y0 < 4 or y1 - prev[i, 1] < 4): s_ += q * (yb - ya) / dt
            flux.append(s_ / L)
            # annihilation: defects of prev with no forward link, +/- pairs within 3 cells
            gone = [i for i in range(len(prev)) if fw[i] < 0]; used = set()
            for i in gone:
                if prev[i, 2] <= 0 or i in used: continue
                for j in gone:
                    if j not in used and prev[j, 2] < 0 and np.hypot((prev[j, 0] - prev[i, 0] + L / 2) % L - L / 2, prev[j, 1] - prev[i, 1]) < 3:
                        ann_y.append(0.5 * (prev[i, 1] + prev[j, 1])); used |= {i, j}; break
            new = [i for i in range(len(f)) if bw[i] < 0]; used = set()
            for i in new:
                if f[i, 2] <= 0 or i in used: continue
                for j in new:
                    if j not in used and f[j, 2] < 0 and np.hypot((f[j, 0] - f[i, 0] + L / 2) % L - L / 2, f[j, 1] - f[i, 1]) < 3:
                        nuc_y.append(0.5 * (f[i, 1] + f[j, 1])); used |= {i, j}; break
        prev = f
    def where(ys):
        ys = np.array(ys); return (np.mean(ys - y0 < 4), np.mean((ys - y0 >= 4) & (y1 - ys >= 4)), np.mean(y1 - ys < 4)) if len(ys) else (np.nan,) * 3
    return dict(run=os.path.basename(run), W=W, Qb=Qbw, Qt=Qtw, act=act, hp=hp / max(nf, 1), hm=hm / max(nf, 1), occ=occ / max(nf, 1),
                Jcross=1000 * (up - down) / (L * T), Jflux=1000 * np.mean(flux) if flux else np.nan, Jflux_se=1000 * np.std(flux) / np.sqrt(max(len(flux), 1)) if flux else np.nan,
                nuc=1000 * len(nuc_y) / (L * T), ann=1000 * len(ann_y) / (L * T), nuc_where=where(nuc_y), ann_where=where(ann_y))
def show(r):
    net = 0.5 * (r["hp"] - r["hm"]); o = r["occ"]
    print(f"{r['run']}: W {r['W']:.0f}, wall charges {r['Qb']:+.0f} (bottom) / {r['Qt']:+.0f} (top), act {r['act']}")
    print(f"   layers: bottom +1/2 {o[0]:.2f}, -1/2 {o[1]:.2f} | top +1/2 {o[2]:.2f}, -1/2 {o[3]:.2f}")
    print("   net charge across the width (bottom -> top): " + " ".join(f"{x:+.2f}" for x in net))
    print(f"   CURRENT (+y, per unit length per 1000 t.u.): mid-line crossings {r['Jcross']:+.3f}; bulk flux {r['Jflux']:+.3f} ± {r['Jflux_se']:.3f}")
    print(f"   nucleation {r['nuc']:.3f}, annihilation {r['ann']:.3f} per unit length per 1000 t.u.; where (bottom band / bulk / top band): nucleation {np.round(r['nuc_where'], 2).tolist()}, annihilation {np.round(r['ann_where'], 2).tolist()}")
if __name__ == "__main__":
    runs = [a for a in sys.argv[1:] if not a.startswith("--")]; res = [analyse(r) for r in runs]
    for r in res: show(r)
    if "--group" in sys.argv:
        G = {}
        for r in res: G.setdefault(re.sub(r"_s\d+$", "", r["run"]), []).append(r)
        print("\nGROUPED: config | n | J_cross | J_flux | nucleation | annihilation")
        for k, v in sorted(G.items()):
            print(f"  {k:40s} | {len(v)} | {np.mean([x['Jcross'] for x in v]):+.3f} | {np.mean([x['Jflux'] for x in v]):+.3f} | {np.mean([x['nuc'] for x in v]):.3f} | {np.mean([x['ann'] for x in v]):.3f}")
