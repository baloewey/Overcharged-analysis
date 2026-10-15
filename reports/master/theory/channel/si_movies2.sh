#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel; F=../../../figstyle/run.sh
while pgrep -f "cavity_movie.py /home/baloewe/Physics/Overcharged/runs/fineflow/eps2.0_s1" >/dev/null; do sleep 20; done
nice $F current_movie.py /home/baloewe/Physics/Overcharged/runs/channel_mixedW/W24_Qb+6_Qt-6_act_s1 /home/baloewe/Physics/Overcharged/runs/channel_narrow/W24_m6_s1 1000 4000 3
nice $F cavity_movie.py /home/baloewe/Physics/Overcharged/runs/swirl_threshold/passive_A0.01_s2 0 600 1 S2_passive_swirl_ejection "Passive circle, k = 6, imposed vortex: the ring of 12 ejects two"
nice $F cavity_movie.py /home/baloewe/Physics/Overcharged/runs/fineflow/eps1.0_s1 1500 2100 2 S1_circle_overcharged "Active circle, k = 6: wall layer of 10, bulk +1, one vortex"
echo done2
