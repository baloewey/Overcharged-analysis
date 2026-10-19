#!/bin/bash
# Wait for the continuum snapshots, then render the tilt-crowding movie (detached; 2026-10-16).
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt
while ps -eo args | grep -q '^python3 active_snapshots'; do sleep 60; done
cat snaps_k4.log snaps_k16.log | tail -4
cd ../tiltcrowd && nice -n 10 ../../../figstyle/run.sh tilt_movie.py 2>&1 | grep -v Warn | tail -3
echo "$(date) movie done"
