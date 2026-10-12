#!/bin/bash
# start the half-integer-k runs once the modulated-winding test (anchmod) has finished; production binary
cd "$(dirname "$0")"
while pgrep -f "mpcd_anchmod.out -i" > /dev/null; do sleep 60; done
EXE=$(realpath ../../mpcd/mpcd.out)
for d in */; do (cd "$d" && setsid -f nohup bash -c "$EXE -i input.json -o . > log.txt 2>&1; grep -q 'Simulation successful' log.txt && touch DONE" < /dev/null > /dev/null 2>&1); sleep 2; done
echo "$(date '+%F %T') started $(ls -d */ | wc -l) runs"
