#!/usr/bin/env python3
"""
Add seeds s04..s10 to the shape sweeps (shapes = +k, shapes_mirror = -k; circle/oval/ellipse; 2026-10-15, user decision:
queue after the mirror sweep, may be cancelled later).

Each new run copies the point's s02/input.json and changes only the seed, continuing the sweeps' own scheme
(seed = base + 100 * i_point + s; base 80000000 for shapes, 81000000 for shapes_mirror; read off s02's seed, so s04..s10
are unique and never collide with s01..s03).
Writes <sweep>/manifest_seeds10.tsv and ~/Overcharged/runs/jobs_shapes_seeds10.txt (+k block first, then -k;
most expensive first within each block).
Usage: extend_shape_seeds.py [--nseeds 10] [--dry-run]
"""
import argparse, glob, json, os
R = os.path.expanduser("~/Overcharged/runs")
ap = argparse.ArgumentParser(); ap.add_argument("--nseeds", type=int, default=10); ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
jobs = []
for sweep in ("shapes", "shapes_mirror"):
    pts = sorted(os.path.dirname(os.path.dirname(p)) for p in glob.glob(os.path.join(R, sweep, "*", "*", "Req*", "k*", "s02", "input.json")))
    rows = []
    for d in pts:
        base = json.load(open(os.path.join(d, "s02", "input.json"))); s0 = base["seed"] - 2
        for s in range(4, a.nseeds + 1):
            rd = os.path.join(d, f"s{s:02d}"); j = dict(base); j["seed"] = s0 + s
            cost = j["species"][0]["pop"] * (j["simSteps"] + j.get("warmUp", 0))
            rows.append(f"{os.path.relpath(rd, os.path.join(R, sweep))}\t{j['seed']}\t{j['species'][0]['pop']}\t{cost}")
            jobs.append((sweep, -cost, rd))
            if not a.dry_run:
                if os.path.exists(os.path.join(rd, "input.json")): raise SystemExit(f"exists, refusing to overwrite: {rd}")
                os.makedirs(rd); json.dump(j, open(os.path.join(rd, "input.json"), "w"), indent=4)
    print(f"{sweep}: {len(pts)} points, {len(rows)} new runs")
    if not a.dry_run:
        open(os.path.join(R, sweep, "manifest_seeds10.tsv"), "w").write("dir\tseed\tpop\tcost\n" + "\n".join(rows) + "\n")
seeds = [json.load(open(p))["seed"] for p in glob.glob(os.path.join(R, "shapes*", "*", "*", "Req*", "k*", "s0[1-3]", "input.json"))]
print("total new runs", len(jobs), "| old seeds", len(seeds), "unique", len(set(seeds)))
if not a.dry_run:
    open(os.path.join(R, "jobs_shapes_seeds10.txt"), "w").write("\n".join(rd for _, _, rd in sorted(jobs)) + "\n")
