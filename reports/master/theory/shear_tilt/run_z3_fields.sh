#!/bin/bash
# zeta 3 rerun with field output (user, 2026-10-19): core depth and wall shear for the d/l_E test; deterministic, so tilts should repeat var_act_3_k*.json
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for k in 4 16; do SAVE=1 TAG=_fields nice -n 10 python3 active_variants.py act $k 3 > var_act_3_k${k}_fields.log 2>&1 & done
wait; echo "$(date) z3 fields done"
