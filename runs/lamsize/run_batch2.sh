#!/bin/bash
# lamsize batch 2 (2026-10-13): 3 more seeds for R_eq 1/1.5/2 x lambda 1.5/1.75/2 and 5 seeds of small circles R_eq 0.5/0.75 x lambda 1.5/1.75/2/3 (67 runs).
# Waits for the after_restart channel chain (bash PID 7052) to finish, then runs 14 in parallel; then summarises per R_eq -> energetics/tiltlaw_lamsize_byR.txt
cd /home/baloewe/Physics/Overcharged/runs/lamsize
# (2026-10-14: started early in parallel with the channel chain, user OK)
echo "$(date) chain finished, starting batch 2"
./run_queue.sh jobs2.txt 14 > queue2.log 2>&1
T=/home/baloewe/Physics/Overcharged/reports/master/theory/energetics; O=$T/tiltlaw_lamsize_byR.txt; : > $O
for r in 0.5 0.75 1.0 1.5 2.0; do echo "R_eq $r" >> $O; OUT=/dev/null nice -n 19 python3 $T/tiltlaw.py Req${r}_*/ 2>/dev/null | tail -n +2 >> $O; done
echo "$(date) batch 2 done and summarised"
