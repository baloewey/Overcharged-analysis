#!/bin/bash
# Talk movies with a fixed frame (user 2026-11-26: S1/S2/S3 jittered because the tight bbox changed with the header text) + new Cassini oval k 6.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel; F=../../../figstyle/run.sh; export MOVIE_OUTDIR=/home/baloewe/Physics/Overcharged/reports/talk_workshop_20261127/movies_fixed
R=/home/baloewe/Physics/Overcharged/runs
nice $F cavity_movie.py $R/fineflow/eps1.0_s1 1500 2100 2 S1_circle_overcharged "Active circle, k = 6: wall layer of 10, bulk +1, one vortex" &
nice $F cavity_movie.py $R/swirl_threshold/passive_A0.01_s2 0 600 1 S2_passive_swirl_ejection "Passive circle, k = 6, imposed vortex: the ring of 12 ejects two" &
nice $F cavity_movie.py $R/fineflow/eps2.0_s1 1500 2100 2 S3_ellipse_leaking "Active ellipse, aspect 2, k = 6: the layer leaks" &
wait; echo done
