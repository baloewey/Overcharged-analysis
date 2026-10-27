#!/bin/bash
# After mirror_esc.py: +k strict escapes, then the circle occupancy-split excursion profiles (item 6).
cd /home/baloewe/Physics/Overcharged/reports/master/theory/crowd_protect/mirror10
while pgrep -f "python3 mirror_esc.py" > /dev/null; do sleep 30; done
NPROC=12 nice -n 19 python3 plus_strict.py > plus_strict.log 2>&1
cd ../drive_barrier && NPROC=12 nice -n 19 python3 circle_occupancy.py > circle_occupancy.log 2>&1
