#!/usr/bin/env python3
"""Provenance (audit 2026-10-19) for two numbers that existed only in RESEARCH.md / a figure legend:
(1) free linear fits tilt = tau0 + B d sqrt(s), active (10-seed circle window: tiltcrowd/defect_flow_main10.txt + wallshear_main10.txt) and passive
    (tiltcrowd/passive_runs.npz), exactly as tilt_report/report_figs.py Fig. 5b; also the depth range by spacing bin (medians).
(2) exit-flag rate per wall defect per 20 t.u. against vacancies (2k-2 minus wall count) from energetics/perdef_ff36.npy, with the run conditions (aspect, lambda).
-> dsqrts_exit.txt"""
import numpy as np
T = "/home/baloewe/Physics/Overcharged/reports/master/theory"
def rows(fn, prefix):
    return [[c.strip() for c in l.split("|")] for l in open(fn) if l.startswith(prefix)]
DFm = rows(f"{T}/tiltcrowd/defect_flow_main10.txt", "rho"); WSm = rows(f"{T}/tiltcrowd/wallshear_main10.txt", "rho"); P = np.load(f"{T}/tiltcrowd/passive_runs.npz")
dd = np.array([float(r[3]) for r in DFm]); ws = {r[0]: float(r[5]) + float(r[6]) for r in WSm}; ds = np.array([ws.get(r[0], np.nan) for r in DFm])
da = np.array([float(r[2]) for r in DFm]); dT = np.array([float(r[6]) for r in DFm]); ok = np.isfinite(dT) & np.isfinite(ds)
L = ["(1) free fits tilt = tau0 + B d sqrt(s) (as tilt report Fig. 5b)"]
for lab, x, y, a, d in (("active", dd[ok] * np.sqrt(ds[ok]), dT[ok], da[ok], dd[ok]), ("passive", P["d"] * np.sqrt(P["s"]), P["tilt"], P["a"], P["d"])):
    B, t0 = np.polyfit(x, y, 1); r = y - (t0 + B * x)
    L.append(f"  {lab}: n {len(x)}; tau0 {t0:.1f} deg, B {B:.1f}; rms residual {r.std():.2f} (raw sd {y.std():.2f})")
    L.append("    depth by spacing bin (median): " + ", ".join(f"[{lo},{hi}) {np.median(d[(a>=lo)&(a<hi)]):.2f}" for lo, hi in ((5, 8), (8, 12), (12, 20), (20, 60)) if ((a>=lo)&(a<hi)).any()))
E = np.load(f"{T}/energetics/perdef_ff36.npy", allow_pickle=True)
L.append("\n(2) exit-flag rate per wall defect per 20 t.u. vs vacancies (perdef_ff36.npy; rows: col 0 frame, col 1 vacancies?, col 2 wall count; see energetics/perdef*.py)")
L.append("  runs: " + ", ".join(sorted({f"eps{e['eps']} lam{e['lam']}" for e in E})) + f"  ({len(E)} runs)")
open(f"/home/baloewe/Physics/Overcharged/reports/master/theory/provenance/dsqrts_exit.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
# --- (2) exit-flag rate vs vacancies, selection as energetics/perdef.py (eps >= 1.3, all lambda in the set) and lambda 2 only
L.append("  exit flag = wall defect crosses beyond 6 cells within 20 t.u. and stays out (perdef.py); vacancies = (2k-2) - wall count in that frame")
for nm, sel in (("eps >= 1.3, all lambda (perdef.py selection)", lambda e: e["eps"] >= 1.3), ("eps >= 1.3, lambda 2", lambda e: e["eps"] >= 1.3 and e["lam"] == 2.0), ("eps 2.0, lambda 2", lambda e: e["eps"] == 2.0 and e["lam"] == 2.0)):
    ok = [e for e in E if "err" not in e and sel(e) and len(e["rows"])]
    R = np.concatenate([np.c_[e["rows"][:, 2], np.full(len(e["rows"]), e["N0"]), e["rows"][:, 6]] for e in ok]); vac = R[:, 1] - R[:, 0]
    L.append(f"  {nm}: {len(ok)} runs, {len(R)} wall-defect frames; rate by vacancies: " + ", ".join(f"{int(v)}: {R[vac == v, 2].mean():.4f} (n {int((vac == v).sum())})" for v in range(0, 7) if (vac == v).sum() > 500))
open(f"/home/baloewe/Physics/Overcharged/reports/master/theory/provenance/dsqrts_exit.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L[-4:]))
