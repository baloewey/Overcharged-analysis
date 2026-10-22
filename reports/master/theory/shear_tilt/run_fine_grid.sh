#!/bin/bash
# Test F (user, 2026-10-22): active one-constant continuum disk on a twice-finer grid (h 0.25; QGRID_H) so that crowded cores (0.3-0.4 cells at h 0.5)
# are resolved; k 4 and 16 at zeta 3 and 7; fields saved (TAG _h025) for the local-field / excess split. Office, nice 19.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
for z in 3 7; do for k in 4 16; do
  QGRID_H=0.25 SAVE=1 TAG=_h025 nice -n 19 python3 active_variants.py act $k $z > var_act_${z}_k${k}_h025.log 2>&1 &
done; done
wait; echo "$(date) fine-grid runs done" >> run_fine_grid.log
