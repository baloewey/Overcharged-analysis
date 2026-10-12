#!/bin/bash
cd "$(dirname "$0")/.."
for act in 0.06 0.11; do
  R="H3high/act$act/s0[123]/k0[2-9] H3high/act$act/s0[123]/k10"
  python3 edge_sites.py $R --no-gr --out analysis_spatial/act${act/./}_sites > analysis_spatial/act${act/./}_sites.txt 2>&1
  python3 circulation.py $R --out analysis_lunch/circ_act${act/./}_3seeds.csv > analysis_lunch/circ_act${act/./}_3seeds.txt 2>&1
done
echo finished
