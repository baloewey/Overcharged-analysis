#!/bin/bash
# Tracker checking movies: 250 frames (500 t.u.) each at 6 fps (~42 s), windows chosen automatically where the linkers disagree.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/crowd_protect/tracker_validation
R=/home/baloewe/Physics/Overcharged/runs
for r in shapes_server/ellipse/eps2.00/Req1.5/k05/s01 shapes_server/ellipse/eps2.00/Req1.5/k16/s01 shapes_server/circle/circle/Req1.5/k13/s01 channel_current/W24_mixed_s3; do
  n=$(echo $r | tr / _); nice -n 19 python3 validation_movie.py $R/$r -1 250 6 > movies/log_$n.txt 2>&1 &
done
wait
