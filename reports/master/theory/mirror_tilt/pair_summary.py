#!/usr/bin/env python3
"""Summary of pair_events.pkl by class (wall density rho bin x R_eq; circles+aspect<1.3 vs aspect>=1.4): per unit wall length and time,
entry rate split arrive/born and exit rate split depart/die/annih; birth-death rates against the excess n: gain(n) = entries at n / time at n,
loss(n) = exits at n / time at n; predicted occupancy ratio P(n+1)/P(n) = gain(n)/loss(n+1) vs measured. Writes pair_summary.txt."""
import pickle, numpy as np
R = [o for o in pickle.load(open("pair_events.pkl", "rb")) if "error" not in o]
bins = [(0.03, 0.05, "0.04"), (0.05, 0.07, "0.06"), (0.08, 0.10, "0.09"), (0.11, 0.15, "0.12-0.14")]
NS = range(-3, 5); out = []
def agg(g):
    occ = {v: sum(o["occ"][v] for o in g) for v in NS}; ev = {}
    for o in g:
        for key, c in o["ev"].items(): ev[key] = ev.get(key, 0) + c
    L = sum(o["P"] * o["T"] for o in g)   # wall length x time
    return occ, ev, L
out.append("Pair events of wall -1/2 (mirror sweep, second half; pair_events.py). Rates per 1000 cell-lengths of wall per 1000 t.u. (entries, exits);")
out.append("gain(n), loss(n): per run per 1000 t.u. at excess n = N - 2|k|. Excess mean = <N> - 2|k|.")
for shape, sel in (("circle+asp<1.3", lambda o: float(o["aspect"]) < 1.3), ("asp>=1.4", lambda o: float(o["aspect"]) >= 1.4)):
    for req in ("1.0", "1.5", "2.0"):
        for lo, hi, nm in bins:
            g = [o for o in R if sel(o) and o["Req_la"].startswith(req) and lo <= float(o["rho"]) < hi]
            if not g: continue
            occ, ev, L = agg(g); tot = lambda typ: sum(ev.get((typ, v), 0) for v in NS)
            ent = tot("arrive") + tot("born"); ex = tot("depart") + tot("die") + tot("annih"); sc = 1e6 / L
            exc = np.mean([o["Nmean"] - 2 * abs(o["k"]) for o in g])
            out.append(f"\n{shape:15s} R_eq {req} rho~{nm:9s} runs {len(g):3d}  excess {exc:5.2f}  entries {ent*sc:6.2f} (arrive {tot('arrive')/max(ent,1):.0%}, born {tot('born')/max(ent,1):.0%})"
                       f"  exits {ex*sc:6.2f} (depart {tot('depart')/max(ex,1):.0%}, annih {tot('annih')/max(ex,1):.0%}, die {tot('die')/max(ex,1):.0%})")
            row = "   n:      " + "".join(f"{v:>8d}" for v in NS)
            row += "\n   P(n):   " + "".join(f"{occ[v]/sum(occ.values()):8.3f}" for v in NS)
            row += "\n   gain:   " + "".join(f"{1000*(ev.get(('arrive',v),0)+ev.get(('born',v),0))/occ[v]*len(g)/len(g):8.2f}" if occ[v] > 200 * len(g) else f"{'-':>8s}" for v in NS)
            row += "\n   loss:   " + "".join(f"{1000*(ev.get(('depart',v),0)+ev.get(('die',v),0)+ev.get(('annih',v),0))/occ[v]:8.2f}" if occ[v] > 200 * len(g) else f"{'-':>8s}" for v in NS)
            out.append(row)
open("pair_summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
