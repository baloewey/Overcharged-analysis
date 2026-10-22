#!/usr/bin/env python3
"""Does the position shift of the wall defects follow the tilt when the tilt changes? (user, 2026-10-22). Mechanism (flow sets the orientation,
anchoring ties orientation to position along the wall): the wall defects sit at the homeotropic points shifted by tilt/(2k-2) in polar angle,
i.e. arg Z_N = tau with Z_N = <exp(i N phi_j)>, N = 2k-2 (anchoring theta = k phi in all runs). Crowded walls tilt less (64 -> 44 deg), and lambda
changes the tilt at fixed crowding, so the predicted shift changes; if positions were fixed by something else the shift would not follow.
Per run: frames with exactly N wall +1/2 (within 4 cells) and coherent orientation (|<e^{i tau}>| > 0.9); circular means of arg Z_N (A) and of
tau (T, signed); D = A - T. Sets: active main sweep (server_main/main_act0.07, k 3-10, all densities), passive imposed vortex (swirl_crowd),
lambda scans (lambda_circle k 6; lambda_crowd k 14). -> lock_shift.txt, lock_shift.npz"""
import glob, json, os, sys, warnings, numpy as np
warnings.filterwarnings('ignore')
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle
R_ = "/home/baloewe/Physics/Overcharged/runs"
SETS = {"active main sweep": sorted(os.path.dirname(p) for p in glob.glob(R_ + "/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/DONE")),
        "passive imposed vortex": sorted(r for r in glob.glob(R_ + "/swirl_crowd/*k*_s*") if os.path.exists(r + "/DONE")),
        "lambda scan k 6": sorted(r for r in glob.glob(R_ + "/lambda_circle/k06_*") if os.path.exists(r + "/DONE")),
        "lambda scan k 14": sorted(r for r in glob.glob(R_ + "/lambda_crowd/k14_*") if os.path.exists(r + "/DONE"))}
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; c = np.array(bc["Q"][:2]); k = bc["customAnch"]; N = int(round(2 * k - 2))
        if N < 2: return None   # k = 1: no wall layer (2k-2 = 0); excluded
        t, fr = read_with_angle(run, c, 0.5); A, T, rr = [], [], []
        for f in fr:
            q = f[f[:, 2] > 0]; ph = np.arctan2(q[:, 1], q[:, 0]); r = np.hypot(q[:, 0], q[:, 1]); m = R - r < 4
            if m.sum() != N: continue
            e = np.exp(1j * (q[m, 3] - ph[m])).mean()
            if abs(e) < 0.9: continue
            A.append(np.exp(1j * np.angle(np.exp(1j * N * ph[m]).mean()))); T.append(np.exp(1j * np.angle(e))); rr.append(r[m].mean())
        if len(A) < 20: return None
        A, T = np.angle(np.mean(A)), np.angle(np.mean(T)); sg = np.sign(T) or 1.0
        return dict(run=os.path.relpath(run, R_), k=k, N=N, spacing=2 * np.pi * np.mean(rr) / N, A=sg * A, T=sg * T, D=sg * np.angle(np.exp(1j * (A - T))), r=np.mean(rr), n=len(rr))
    except Exception:
        return None
if __name__ == "__main__":
    L, OUT = ["set | runs | slope of A on T (rad/rad) | median D (rad) | |D| 90th pct | median shift error in cells (D/N * r) | T range (deg)"], {}
    for nm, runs in SETS.items():
        with Pool(8) as p: allres = p.map(one, runs, chunksize=2)
        res = [x for x in allres if x]; L.append(f"  ({nm}: {len(runs)} runs, {len(runs) - len(res)} excluded: k = 1 or fewer than 20 coherent frames with exactly 2k-2 wall +1/2)")
        if not res: L.append(f"{nm}: no runs"); continue
        A = np.array([x["A"] for x in res]); T = np.array([x["T"] for x in res]); D = np.array([x["D"] for x in res]); sp = np.array([x["spacing"] for x in res])
        err = np.array([x["D"] / x["N"] * x["r"] for x in res]); sl = np.polyfit(T, A, 1)[0] if np.ptp(T) > 0.05 else np.nan
        OUT[nm] = dict(A=A, T=T, D=D, sp=sp, k=np.array([x["k"] for x in res]))
        L.append(f"{nm} | {len(res)} | {sl:.2f} | {np.median(D):+.3f} | {np.percentile(np.abs(D), 90):.3f} | {np.median(err):+.3f} | {np.degrees(T.min()):.1f}-{np.degrees(T.max()):.1f}")
        for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
            m = (sp >= lo) & (sp < hi)
            if m.sum(): L.append(f"   spacing [{lo},{hi}): {m.sum():3d} runs, tilt {np.degrees(np.median(T[m])):5.1f} deg, shift x N = arg Z {np.degrees(np.median(A[m])):5.1f} deg, D median {np.degrees(np.median(D[m])):+5.1f} deg")
    np.savez("lock_shift.npz", **{f"{nm}|{q}": v for nm, d in OUT.items() for q, v in d.items()})
    open("lock_shift.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
