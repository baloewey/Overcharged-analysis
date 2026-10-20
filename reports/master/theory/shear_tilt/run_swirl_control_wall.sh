#!/bin/bash
# shear-only control, batch 2 (2026-10-20): body force EVERYWHERE in the disk (SWIRL_FRAC 1.05; batch 1 forced only r < 0.6R and raised the wall shear only 0.55 -> 0.7).
# Target the zeta-8 wall shear (~1.9-2.1). zeta 3 (laminar) + A 0.02/0.05/0.1/0.2, k 4/16, fields saved, outputs tagged _wall.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for A in 0.02 0.05 0.1 0.2; do for k in 4 16; do SWIRL_FRAC=1.05 SAVE=1 TAG=_wall nice -n 10 python3 active_variants.py swirl $k $A > var_swirl_${A}_k${k}_wall.log 2>&1 & done; done
wait; echo "$(date) swirl control (wall) done"
