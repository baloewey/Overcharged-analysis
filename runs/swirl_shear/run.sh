#!/bin/bash
# Passive imposed-vortex crowding test (2026-10-15): 12 runs at once, low priority; ~2-3 h each.
cd "$(dirname "$0")"; B=/home/baloewe/Physics/Overcharged/bin/mpcd_swirl.out
for d in ${PAT:-k*_A*}/; do (cd "$d" && setsid -f nohup nice -n 5 bash -c "$B -i input.json -o . > log.txt 2>&1; grep -q 'Simulation successful' log.txt && gzip -f flowfield.dat defects.dat; touch DONE" < /dev/null > /dev/null 2>&1); sleep 2; done
echo "$(date '+%F %T') started runs ${PAT:-k*_A*}"
