"""Effective restoring force on wall +1/2 vs depth (from selfprop_events.npz): per depth bin, drift <dd> and variance <dd^2> - <dd>^2 of the depth change over 10 t.u.
(dd > 0 = moving away from the wall). Fokker-Planck (Ito): U'(d) / kT_eff = -2 <dd> / var(dd)  ->  U(d) by integration, relative to its minimum. Also the 'noise' var(dd).
Classes: circle-like vs elongated, by R_eq; elongated by region. Writes restoring.txt."""
import numpy as np
Z = np.load("selfprop_events.npz"); D, DR, EL, RQ, RG = Z["depth"], Z["drift"], Z["elong"], Z["Req"], Z["region"]; ok = np.isfinite(DR)
edges = np.arange(0.0, 4.01, 0.5); c = 0.5 * (edges[1:] + edges[:-1]); L = []
def prof(m):
    mu, var, n = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = m & ok & (D >= lo) & (D < hi); x = DR[k]; n.append(k.sum()); mu.append(x.mean() if len(x) > 100 else np.nan); var.append(x.var() if len(x) > 100 else np.nan)
    mu, var = np.array(mu), np.array(var); F = -2 * mu / var; U = np.concatenate([[0], np.cumsum(F[:-1] * 0.5 * 0 + 0.5 * (F[:-1] + F[1:]) * 0.5)])  # trapezoid over bin centres
    return mu, var, F, U - np.nanmin(U), np.array(n)
groups = [(f"{'elongated' if e else 'circle-like'} R_eq {R}", (EL == e) & (RQ == R)) for e in (False, True) for R in (1.0, 1.5, 2.0)] + \
         [(f"elongated {['tip', 'middle', 'flank'][g]}", EL & (RG == g)) for g in (0, 1, 2)]
L.append("depth bins (centre): " + " ".join(f"{x:5.2f}" for x in c))
for name, m in groups:
    mu, var, F, U, n = prof(m)
    L.append(f"\n{name}  (n {m.sum()})")
    L.append("   drift <dd> (cells/10 t.u.): " + " ".join(f"{x:+6.3f}" for x in mu))
    L.append("   variance  (cells^2)       : " + " ".join(f"{x:6.3f}" for x in var))
    L.append("   force -U'/kT (1/cell)     : " + " ".join(f"{-x:+6.2f}" for x in F))
    L.append("   U(d)/kT_eff               : " + " ".join(f"{x:6.2f}" for x in U))
print("\n".join(L)); open("restoring.txt", "w").write("\n".join(L) + "\n")
