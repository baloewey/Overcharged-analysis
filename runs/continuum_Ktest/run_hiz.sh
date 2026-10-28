#!/bin/bash
# Continuum K test at higher activity (2026-10-28, user: does the crowded continuum layer leak once the bulk is turbulent?).
# Same geometry as run.sh (aspect-2 ellipse, R_eq 31.83, lambda 2), k 5 / 14. dt as the circle zeta scan (active_variants.py): 0.005 at zeta 3, 0.002 at zeta 8.
#   z3   zeta 3, noise 0.005, 3 seeds          (laminar-ish, crowding drop appears in circles)
#   z8   zeta 8, noise 0.005, 3 seeds          (turbulent bulk with pair creation, vortex kept in circles)
#   z8N05 zeta 8, noise /sqrt2, 2 seeds        (thermal-activation / matched-l_a K test at high activity)
cd "$(dirname "$0")"; PY=$HOME/Overcharged/venv_analysis/bin/python; S=$HOME/Overcharged/continuum/ellipse_active_K.py
launch() { o=$1; shift
  setsid -f nohup nice -n 19 bash -c "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY $S --eps 2 --k $K --T0 600 --T1 4000 --seed $SD --track 2 $* --out $o.json > $o.log 2>&1; date > $o.DONE" < /dev/null > /dev/null 2>&1; }
for K in 5 14; do
  for SD in 1 2 3; do launch z3_k${K}_s${SD} --zeta 3 --dt 0.005; launch z8_k${K}_s${SD} --zeta 8 --dt 0.002; done
  for SD in 1 2; do launch z8N05_k${K}_s${SD} --zeta 8 --dt 0.002 --noise 0.0035355; done
done
echo "$(date '+%F %T') started 16 high-activity runs"
