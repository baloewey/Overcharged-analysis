"""Second activity scan (first one: zeta <= 0.03 gave no flow). zeta 0.1-1 at lam 2 and 1.2; lam 4 with dt 0.01; passive+bf stronger body force (ejection?).
-> scan_active2.json / .txt (same format as scan_active.py)."""
import json, numpy as np
from multiprocessing import Pool
import active_bf as B
np.seterr(all="ignore")
JOBS = [("active", 2.0, z, 0.0, 0.02) for z in (0.1, 0.3, 1.0)] + [("active", 1.2, z, 0.0, 0.02) for z in (0.1, 0.3)] + [("active", 4.0, 0.3, 0.0, 0.005)] + \
       [("passive+bf", 2.0, 0.0, A, 0.02) for A in (0.03, 0.1)] + [("passive+bf", 4.0, 0.0, 0.03, 0.005)]
def one(j):
    tag, lam, zeta, A, dt = j; obs, (q1, q2, vx, vy) = B.run(lam, zeta, A, ring=12, T0=300, T1=4000, dt=dt, every=250, nflow=1); last = obs[-6:]
    f = lambda key: float(np.nanmean([o.get(key, np.nan) for o in last]))
    zs = [np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]; z = np.nanmean(zs) if zs else np.nan
    np.savez_compressed(f"state2_{tag}_lam{lam}_z{zeta}_A{A}.npz", q1=q1, q2=q2, vx=vx, vy=vy)
    return dict(tag=tag, lam=lam, zeta=zeta, A=A, dt=dt, nwall=f("nwall"), nbulk_p=f("nbulk_p"), nbulk_m=f("nbulk_m"), coh=f("coh"),
                tilt=float(np.degrees(np.angle(z))) if zs else np.nan, absc=float(np.nanmean([abs(o["circ"]) for o in last])), vrms=f("vrms"), obs=obs)
if __name__ == "__main__":
    with Pool(9) as p: res = p.map(one, JOBS, chunksize=1)
    json.dump(res, open("scan_active2.json", "w"), default=float); L = []
    for r in res:
        L.append(f"{r['tag']:10s} lam {r['lam']:3.1f} zeta {r['zeta']:.2f} A {r['A']:.3f} dt {r['dt']}: wall {r['nwall']:4.1f} bulk +{r['nbulk_p']:.1f}/-{r['nbulk_m']:.1f} "
                 f"coh {r['coh']:.2f} tilt {r['tilt']:+6.1f} |circ| {r['absc']:.2f} vrms {r['vrms']:.1e}")
    open("scan_active2.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
