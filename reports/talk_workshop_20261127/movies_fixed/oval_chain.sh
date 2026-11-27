#!/bin/bash
# Render the narrow-neck oval movies once runs/oval_neck s1 runs finish (user 2026-11-26).
R=/home/baloewe/Physics/Overcharged/runs/oval_neck
until [ -f $R/ba1.10_k06_s1/DONE ] && [ -f $R/ba1.05_k06_s1/DONE ]; do sleep 60; done
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel; F=../../../figstyle/run.sh; export MOVIE_OUTDIR=/home/baloewe/Physics/Overcharged/reports/talk_workshop_20261127/movies_fixed
nice $F cavity_movie.py $R/ba1.10_k06_s1 1600 2200 2 oval_neck_ba1.10_k6 "Active Cassini oval with a neck (b/a 1.10, aspect 2.5), k = 6" &
nice $F cavity_movie.py $R/ba1.05_k06_s1 1600 2200 2 oval_neck_ba1.05_k6 "Active Cassini oval, narrow neck (b/a 1.05, aspect 2.6), k = 6" &
wait; echo done
