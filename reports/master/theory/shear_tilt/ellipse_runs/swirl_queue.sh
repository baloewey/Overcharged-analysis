#!/bin/bash
# Continuum forced-vortex + localized-activity set (2026-10-10): launch each job when fewer than 12 continuum runs are active.
J=(); for e in 1.0 1.6 2.0; do
  J+=("sw_passive_eps$e|--eps $e --zeta 0 --swirl 0.1")
  J+=("sw_band_eps$e|--eps $e --zeta 1.0 --bandonly 4 --swirl 0.1")
  J+=("sw_core_eps$e|--eps $e --zeta 1.0 --actmode core --swirl 0.1")
  J+=("sw_bulk_eps$e|--eps $e --zeta 1.0 --nowall 4 --swirl 0.1")
  J+=("sw_all_eps$e|--eps $e --zeta 1.0 --swirl 0.1")
done
for j in "${J[@]}"; do n=${j%%|*}; a=${j#*|}
  while [ $(ps -eo args | grep -c "^python3 ../ellipse_active.py") -ge 12 ]; do sleep 120; done
  echo "$(date) start $n" >> swirl_queue.log
  setsid -f nohup bash -c "python3 ../ellipse_active.py $a --T1 10000 --out $n.json > $n.log 2>&1"; sleep 5
done
echo "$(date) all launched" >> swirl_queue.log
