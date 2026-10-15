#!/bin/bash
# Start the 10-seed circle dynamics after the sweep-table rerun (runs/cavity/rerun_seeds10.sh, pid 169739) exits (2026-10-15).
cd /home/baloewe/Physics/Overcharged/reports/master/theory/circles10
while kill -0 169739 2>/dev/null; do sleep 120; done
echo "$(date) sweep tables done; starting circles10_dynamics"
nice -n 19 python3 circles10_dynamics.py 10
echo "$(date) circles10_dynamics exited $?"
