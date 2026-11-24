#!/bin/bash
# Continuum band/bulk tilt with and without MPCD's chi (2026-11-23, user): circle R 31.83 (eps 1), lambda 2, noise 0.005, k 6 / 14, chi 1 / 0.5,
# zeta 1 (laminar single vortex; dt 0.02, 3000 t.u.) and zeta 8 (turbulent bulk with pair creation, vortex kept; dt 0.002, 2000 t.u.), 2 seeds.
# --profile 1: ring profiles of the director angle from the tangent and of the azimuthal velocity vs depth every 100 t.u.
cd "$(dirname "$0")"; PY=$HOME/Overcharged/venv_analysis/bin/python; S=$HOME/Overcharged/continuum/ellipse_active_K.py
for z in "1|0.02|3000" "8|0.002|2000"; do IFS="|" read zeta dt T1 <<< "$z"
for k in 6 14; do for chi in 1 0.5; do for s in 1 2; do
  o=snap_z${zeta}_k${k}_chi${chi}_s${s}
  setsid -f nohup nice -n 19 bash -c "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY $S --eps 1 --zeta $zeta --k $k --T0 600 --T1 $T1 --dt $dt --seed $s --chi $chi --profile 1 --snap 1 --out $o.json > $o.log 2>&1; date > $o.DONE" < /dev/null > /dev/null 2>&1
done; done; done; done
echo "$(date '+%F %T') started 16 continuum chi runs with snapshots"
