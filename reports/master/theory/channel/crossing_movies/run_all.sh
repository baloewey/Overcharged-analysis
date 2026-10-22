#!/bin/bash
# Render the five crossing movies in parallel (nice 19).
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel
for r in channel_mixedW/W24_Qb+6_Qt-6_act_s1 channel_mixedW/W60_Qb+6_Qt-6_act_s1 channel_current/W24_mixed_s3 channel_current/W24_sym_s3 channel_fix/crowd_W40_m18_mixed_s1; do
  n=$(basename $r); nice -n 19 python3 crossing_movie.py /home/baloewe/Physics/Overcharged/runs/$r 40 25 > crossing_movies/log_$n.txt 2>&1 &
done
wait
