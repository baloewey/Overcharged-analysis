#!/bin/bash
# Narrow wall-band knockout (d = 2.5) after the swirl queue has launched everything.
while pgrep -f "swirl_queue.sh" > /dev/null; do sleep 120; done
for e in 1.0 1.2 1.4 1.6 2.0 2.5; do n=narrow_eps${e}_z1.0
  while [ $(ps -eo args | grep -c "^python3 ../ellipse_active.py") -ge 12 ]; do sleep 120; done
  echo "$(date) start $n" >> swirl_queue.log
  setsid -f nohup bash -c "python3 ../ellipse_active.py --eps $e --zeta 1.0 --nowall 2.5 --T1 12000 --out $n.json > $n.log 2>&1"; sleep 5
done
