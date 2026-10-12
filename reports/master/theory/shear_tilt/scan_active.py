"""Activity + backflow scan (active_bf.py) -> scan_active.json / scan_active.txt. Steady = last 8 samples (2000 t.u.)."""
import json, numpy as np
from multiprocessing import Pool
import active_bf as B
np.seterr(all="ignore")
JOBS = [("active", 2.0, z, 0.0) for z in (0.002, 0.004, 0.01, 0.03)] + [("active", l, z, 0.0) for z in (0.004, 0.01) for l in (0.5, 1.2, 4.0)] + \
       [("passive+bf", l, 0.0, A) for l in (2.0, 4.0) for A in (0.003, 0.01)]
def one(j):
    tag, lam, zeta, A = j; obs, (q1, q2, vx, vy) = B.run(lam, zeta, A, ring=12, T0=300, T1=6000, dt=0.04, every=250, nflow=1); last = obs[-8:]
    f = lambda key: float(np.nanmean([o.get(key, np.nan) for o in last]))
    z = np.nanmean([np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]) if any("tilt" in o for o in last) else np.nan
    np.savez_compressed(f"state_{tag}_lam{lam}_z{zeta}_A{A}.npz", q1=q1, q2=q2, vx=vx, vy=vy)
    return dict(tag=tag, lam=lam, zeta=zeta, A=A, nwall=f("nwall"), nbulk_p=f("nbulk_p"), nbulk_m=f("nbulk_m"), coh=f("coh"), tilt=float(np.degrees(np.angle(z))) if np.isfinite(abs(z)) else np.nan,
                circ=f("circ"), absc=float(np.nanmean([abs(o["circ"]) for o in last])), vrms=f("vrms"), obs=obs)
if __name__ == "__main__":
    with Pool(7) as p: res = p.map(one, JOBS, chunksize=1)
    json.dump(res, open("scan_active.json", "w"), default=float); L = []
    for r in res:
        L.append(f"{r['tag']:10s} lam {r['lam']:3.1f} zeta {r['zeta']:.3f} (R/l_a {31.83 / np.sqrt(1.535 / r['zeta']) if r['zeta'] else 0:.2f}) A {r['A']:.3f}: wall {r['nwall']:4.1f} "
                 f"bulk +{r['nbulk_p']:.1f}/-{r['nbulk_m']:.1f} coh {r['coh']:.2f} tilt {r['tilt']:+6.1f} |circ| {r['absc']:.2f} vrms {r['vrms']:.1e}")
    open("scan_active.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
