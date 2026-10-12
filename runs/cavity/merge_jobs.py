#!/usr/bin/env python3
"""Merge manifests into one job list, most expensive first.  merge_jobs.py OUT.txt manifest1.csv [...]"""
import csv, os, sys
out, mans = sys.argv[1], sys.argv[2:]
rows = [r for m in mans for r in csv.DictReader(open(m))]
rows.sort(key=lambda r: -float(r["cpu_h"]))
with open(out, "w") as f:
    f.writelines(os.path.abspath(r["dir"]) + "\n" for r in rows)
ch = sum(float(r["cpu_h"]) for r in rows)
gb = sum(float(r["storage_MB"]) for r in rows) / 1e3
big = max(float(r["cpu_h"]) for r in rows)
print(f"{out}: {len(rows)} jobs, {ch:.0f} office-single core-h (~{ch/22:.0f} h on the server at 96-way), "
      f"fields {gb:.0f} GB uncompressed, largest {big:.1f} core-h (~{big*4.4:.0f} h at 96-way)")
