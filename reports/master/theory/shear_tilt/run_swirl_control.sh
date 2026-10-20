#!/bin/bash
# shear-only control (user, 2026-10-20): laminar zeta 3 disk + imposed azimuthal body force A; does the crowding drop reach the zeta-8 value (10.7 deg) at the
# zeta-8 wall shear (1.9-2.1) without turbulence? k 4 and 16, A 0.005/0.01/0.02/0.04, fields saved. 8 at once, detached.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for A in 0.005 0.01 0.02 0.04; do for k in 4 16; do SAVE=1 nice -n 10 python3 active_variants.py swirl $k $A > var_swirl_${A}_k$k.log 2>&1 & done; done
wait; echo "$(date) swirl control done"
