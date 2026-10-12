"""Summary of the active continuum ellipse runs (second half of each run): wall +1/2, bulk charge, tilt coherence, tilt signed with the circulation, |circulation|, P(two vortices).
Compare with the MPCD shape sweep (phase_shapes): transition at aspect 1.25-1.37, coherence ~0.45 and two vortices 32-84% above it, vacancies 1-4.5."""
import glob, json, numpy as np
rows = []
for f in sorted(glob.glob("eps*_z*.json")):
    J = json.load(open(f)); a = J["args"]; o = J["obs"]; o = o[len(o) // 2:]
    g = lambda k: np.nanmean([x.get(k, np.nan) for x in o])
    st = np.nanmean([x.get("tilt", np.nan) * np.sign(x["circ"]) for x in o if "tilt" in x])
    rows.append((a["eps"], a["zeta"], g("nwall"), 0.5 * (g("nbp") - g("nbm")), g("coh"), st, np.nanmean([abs(x["circ"]) for x in o]), g("two")))
print(f"{'eps':>4s} {'zeta':>5s} {'wall':>5s} {'vac':>5s} {'Qbulk':>6s} {'coh':>5s} {'tilt':>6s} {'|circ|':>6s} {'P(two)':>6s}")
for e, z, w, q, c, t, ci, two in sorted(rows):
    print(f"{e:4.1f} {z:5.1f} {w:5.2f} {10 - w:5.2f} {q:+6.2f} {c:5.2f} {t:+6.1f} {ci:6.2f} {two:6.2f}")
