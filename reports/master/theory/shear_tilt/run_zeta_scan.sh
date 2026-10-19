#!/bin/bash
# continuum active disk, intermediate activities (user, 2026-10-19): want bulk turbulence (pair creation, disordered spiral) while the crowded k 16 disk keeps its
# vortex and shows a crowding drop. zeta 4-8, k 4 and 16, fields saved (d/l_E). 6 at a time, detached.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
jobs_list=""; for z in 4 5 6 7 8; do for k in 4 16; do jobs_list="$jobs_list $z:$k"; done; done
printf "%s\n" $jobs_list | xargs -P 6 -I{} bash -c 'z=${0%%:*}; k=${0##*:}; SAVE=1 nice -n 10 python3 active_variants.py act $k $z > var_act_${z}_k${k}.log 2>&1' {}
echo "$(date) zeta scan done"
