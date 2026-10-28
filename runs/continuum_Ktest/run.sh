#!/bin/bash
# Continuum K test (2026-10-27, user): aspect-2 ellipse (R_eq 31.83, lambda 2, zeta 1, noise 0.005), k 5 / 14, 2 seeds, 600 + 4000 t.u., defects every 2 t.u.
#   ref    Lc 1, G 1                      reference
#   K2     Lc 2                           K alone doubled (l_a, core change)
#   K2G05  Lc 2, G 0.5                    K and gamma_1 doubled together (as raising mfpot in MPCD: rot_visc.txt)
#   N05    noise 0.005/sqrt2              = K, zeta, a, c all doubled at fixed eta, Gamma up to a time rescale: only kBT/K halves (thermal activation test)
cd "$(dirname "$0")"; PY=$HOME/Overcharged/venv_analysis/bin/python; S=$HOME/Overcharged/continuum/ellipse_active_K.py
for k in 5 14; do for s in 1 2; do
  for c in "ref|" "K2|--Lc 2" "K2G05|--Lc 2 --G 0.5" "N05|--noise 0.0035355"; do
    tag=${c%%|*}; extra=${c#*|}; o=${tag}_k${k}_s${s}
    setsid -f nohup nice -n 19 bash -c "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY $S --eps 2 --zeta 1 --k $k --T0 600 --T1 4000 --seed $s --track 2 $extra --out $o.json > $o.log 2>&1; date > $o.DONE" < /dev/null > /dev/null 2>&1
  done; done; done
echo "$(date '+%F %T') started 16 runs"
