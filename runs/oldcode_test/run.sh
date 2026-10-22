#!/bin/bash
# Activity-scale shift test (2026-10-22): Jan-2024 code (be08195, built in ~/Overcharged/oldcode/mpcd_be08195) vs production, same inputs; nice 19.
D=/home/baloewe/Overcharged/runs/oldcode_test; cd $D
declare -A EXE=( [old]=/home/baloewe/Overcharged/oldcode/mpcd_be08195/mpcd.out [new]=/home/baloewe/Overcharged/mpcd/mpcd.out )
while read code name; do
  ( cd $name; echo "$(date '+%F %T') start $name" >> ../queue.log
    nice -n 19 ${EXE[$code]} -i input.json -o . > log.txt 2>&1 & echo $! > pid; wait $!
    grep -q "Simulation successful" log.txt && { gzip -f defects.dat directorfield.dat log.txt 2>/dev/null; touch DONE; }
    echo "$(date '+%F %T') done  $name" >> ../queue.log ) &
done < jobs.txt
wait; echo "$(date '+%F %T') all finished" >> queue.log
