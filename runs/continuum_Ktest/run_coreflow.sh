#!/bin/bash
# Continuum flow-at-core test (2026-10-28, user): as ref (zeta 1, dt 0.02) and z8 (zeta 8, dt 0.002), aspect-2 ellipse k 5 / 14, 2 seeds, with --coreflow 1
# (defect rows x, y, q, vx, vy; flow = disk mean r < 2 around the core). MPCD analogue: runs/coreflow.
cd "$(dirname "$0")"; PY=$HOME/Overcharged/venv_analysis/bin/python; S=$HOME/Overcharged/continuum/ellipse_active_K.py
for k in 5 14; do for s in 1 2; do for c in "cf1|--zeta 1 --dt 0.02" "cf8|--zeta 8 --dt 0.002"; do
  tag=${c%%|*}; extra=${c#*|}; o=${tag}_k${k}_s${s}
  setsid -f nohup nice -n 19 bash -c "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY $S --eps 2 --k $k --T0 600 --T1 4000 --seed $s --track 2 --coreflow 1 $extra --out $o.json > $o.log 2>&1; date > $o.DONE" < /dev/null > /dev/null 2>&1
done; done; done
echo "$(date '+%F %T') started 8 coreflow runs"
