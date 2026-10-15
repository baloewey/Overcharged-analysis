#!/bin/bash
# Candidate SI movies (2026-10-15). Runs sequentially, nice.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel
F=../../../figstyle/run.sh
nice $F cavity_movie.py /home/baloewe/Physics/Overcharged/runs/fineflow/eps1.0_s1 1500 2100 2 S1_circle_overcharged "Active circle, k = 6: wall layer of 10, bulk +1, one vortex"
nice $F cavity_movie.py /home/baloewe/Physics/Overcharged/runs/swirl_threshold/passive_A0.01_s2 0 600 1 S2_passive_swirl_ejection "Passive circle, k = 6, imposed vortex: the ring of 12 ejects two"
nice $F cavity_movie.py /home/baloewe/Physics/Overcharged/runs/fineflow/eps2.0_s1 1500 2100 2 S3_ellipse_leaking "Active ellipse, aspect 2, k = 6: the layer leaks"
echo done
