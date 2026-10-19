#!/bin/bash
# continuum noise/activity variants (2026-10-19), detached
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for k in 4 16; do
  for D in 0.002 0.01; do nice -n 10 python3 active_variants.py noise $k $D > var_noise_${D}_k$k.log 2>&1 & done
  for Z in 3 10; do nice -n 10 python3 active_variants.py act $k $Z > var_act_${Z}_k$k.log 2>&1 & done
done
wait; echo "$(date) variants done"
