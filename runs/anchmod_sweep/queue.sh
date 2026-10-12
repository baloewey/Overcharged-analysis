#!/bin/bash
# start the site-arrangement sweep after the modulated-winding test and the half-integer runs have finished (branch binary)
cd "$(dirname "$0")"
sleep 300
while ps -eo args | grep -q "^/home/baloewe/Physics/Overcharged/bin/mpcd_anchmod.out" || [ ! -s ../halfk/queue.log ] || ps -eo args | grep -q "^/home/baloewe/Physics/Overcharged/mpcd/mpcd.out"; do sleep 120; done
B=/home/baloewe/Physics/Overcharged/bin/mpcd_anchmod.out
for d in */; do (cd "$d" && setsid -f nohup bash -c "$B -i input.json -o . > log.txt 2>&1; grep -q 'Simulation successful' log.txt && touch DONE" < /dev/null > /dev/null 2>&1); sleep 2; done
echo "$(date '+%F %T') started $(ls -d */ | wc -l) runs"
