#!/bin/bash
# Continuum core-size test (2026-10-28, user, option 1): change the core xi = sqrt(Lc/a) at fixed K and zeta (fixed l_a), S0 kept (c = 2a/S0^2).
# Same geometry and settings as run.sh ref (aspect-2 ellipse, zeta 1, noise 0.005, dt 0.02), k 5 / 14, 2 seeds.
#   a05  a 0.5 (xi 1.41, core x sqrt2)      a2  a 2 (xi 0.71 = 1.4 grid cells, near resolution)
cd "$(dirname "$0")"; PY=$HOME/Overcharged/venv_analysis/bin/python; S=$HOME/Overcharged/continuum/ellipse_active_K.py
for k in 5 14; do for s in 1 2; do for A in 0.5 2; do
  o=a${A/./}_k${k}_s${s}
  setsid -f nohup nice -n 19 bash -c "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY $S --eps 2 --zeta 1 --k $k --T0 600 --T1 4000 --seed $s --track 2 --a $A --out $o.json > $o.log 2>&1; date > $o.DONE" < /dev/null > /dev/null 2>&1
done; done; done
echo "$(date '+%F %T') started 8 core-size runs"
