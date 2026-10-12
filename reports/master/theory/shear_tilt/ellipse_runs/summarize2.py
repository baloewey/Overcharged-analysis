"""Continuum ellipses, second half: base (eps*_z*), wall-band activity off (nowall_*), mirror k = -6 (mirror_*). Vacancies relative to 2|k| - 2 = 10 (k = 6) or 2|k| + 2 = 14 (k = -6)."""
import glob, json, numpy as np
def summ(f):
    J = json.load(open(f)); a = J["args"]; o = J["obs"][len(J["obs"]) // 2:]; k = a.get("k", 6); full = 2 * abs(k) - 2 if k > 0 else 2 * abs(k) + 2
    w = np.mean([x["nwall"] for x in o]); qb = 0.5 * (np.mean([x["nbp"] for x in o]) - np.mean([x["nbm"] for x in o]))
    v = np.mean([x["vrms"] for x in o])
    return a["eps"], a["zeta"], w, full - w, qb, np.nanmean([x.get("coh", np.nan) for x in o]), np.mean([abs(x["circ"]) for x in o]), np.mean([x["two"] for x in o]), v
for lab, pat in (("base (k 6)", "eps*_z1.0.json"), ("big R_eq 47.7 (k 9)", "big_eps*_z1.0.json"), ("wall band inactive (k 6)", "nowall_eps*_z1.0.json"), ("mirror (k -6)", "mirror_eps*_z1.0.json"), ("matched flow, uniform zeta 0.75", "match_eps*_z0.75.json"), ("narrow band off (2.5)", "narrow_eps*_z1.0.json"), ("swirl 0.1, passive", "sw_passive_eps*.json"), ("swirl, wall band only active", "sw_band_eps*.json"), ("swirl, cores only active", "sw_core_eps*.json"), ("swirl, bulk only active", "sw_bulk_eps*.json"), ("swirl, active everywhere", "sw_all_eps*.json")):
    print(f"== {lab}:  eps  wall  vac  Qbulk  coh  |circ|  P(two)  v_rms")
    for r in sorted(summ(f) for f in glob.glob(pat)):
        print(f"   {r[0]:4.1f}  {r[2]:5.2f} {r[3]:5.2f} {r[4]:+5.2f} {r[5]:5.2f} {r[6]:5.2f} {r[7]:5.2f} {r[8]:6.3f}")
