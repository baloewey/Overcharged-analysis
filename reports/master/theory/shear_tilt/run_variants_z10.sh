#!/bin/bash
# rerun of the two zeta-10 continuum variants lost in the power cut (2026-10-19 ~13:10)
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for k in 4 16; do nice -n 10 python3 active_variants.py act $k 10 > var_act_10_k$k.log 2>&1 & done
wait; echo "$(date) z10 variants done"
