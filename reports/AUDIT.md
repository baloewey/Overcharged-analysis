# Periodic audit procedure (started 2026-10-15)

User decision: audit our progress periodically. Log every audit in RESEARCH.md ("AUDIT <date>") with the list of changes.

## When
- Before any draft goes to Tyler or out of the group.
- After any new data batch that feeds a quoted number (e.g. 10-seed reruns, mirror sweep, shape 10 seeds).
- At least once a week otherwise.

## What
1. **Numbers.** Every number in `paper/prl/main.tex`, `si.tex`, `paper/long/main.tex`, `theory_report/theory_report.tex` traced to its primary output file (not to another document). Check:
   - value and range; ranges must not mix parameter values (k, shape, size) while being stated for one;
   - stated conditions (k, R/ℓa, λ, ρ, aspect, seeds) match the source;
   - units and definitions (per length vs per area; coherence vs mean angle; kBT; t.u.);
   - the same number agrees across the four documents.
   Keep `paper/prl/NUMBERS.md` (one row per number: source file, conditions) current; numbers without an output file get a summary file in `master/theory/provenance/`.
2. **Columns and definitions.** Any table column a claim relies on is checked against the script that wrote it (lesson: `flow_tilt.csv` `tilt_coh` is the mean angle, not the coherence).
3. **Claims vs data scope.** Does a general statement rest on a single condition (lesson: "≈62°" was one sparse-layer case; the tilt falls in crowded layers)?
4. **Superseded results.** Retractions and corrections in RESEARCH.md propagated to every document.
5. **Figures.** Captions match panels and data; figures regenerated after data changes.
6. **Open TODOs.** Red TODOs in the drafts listed with status.

## How
A read-only agent pass (as on 2026-10-15) producing a table: claim | source | value/conditions | verdict (OK / MISMATCH / UNSOURCED / AMBIGUOUS) | fix. Fixes applied after the user's review; NUMBERS.md and RESEARCH.md updated; commit.

## Log
- 2026-10-15: first PRL audit (6 mismatches, 2 unsourced, ~10 ambiguous; all fixed). Second pass (10-seed circle numbers) pending.
- 2026-10-16: second pass (10-seed circle numbers): all unchanged within s.e.; 3 fixes applied (k range wording, spacing range, mirror band definition).
- 2026-10-19: audit of everything added that day (theory report, PRL v1/v2, SI, tilt report; agent pass → reports/audit_20261019.md): 4 mismatches (SI −½ gain ranges 1.98–2.06 → 1.8–2.1 and 0.1–0.2 → 0.1–0.3; '−½ tilt a third' → about half; ζ 3 continuum sentence out of date), 2 unsourced (active d√s fit, exit rates → provenance/dsqrts_exit.py), 13 ambiguous (scope/wording: −½ 'adds two' → 'up to two', traffic restricted to R_eq ≥ 1.5, annihilation vs lost tracks, story-box drop range, 'within 3°' → '≈3°', restoring force ellipses only, depths, ring lean, turbulent wording, activity robustness scope, high-λ claim scope). All applied after the user's OK; backups *.tex.preaudit.
- 2026-10-20: audit of the 20 Oct additions (continuum ζ scan, density 40, seeded mirror test; reports/audit_20261020.md): numbers all correct; 3 interpretation mismatches + 10 ambiguous wordings, all applied after the user's OK.
