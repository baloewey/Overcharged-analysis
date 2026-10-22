#!/bin/bash
# Kcalib shear-susceptibility runs (2026-10-22): 4 runs, all at once, nice 19 (server shared with the shape queue). DONE marker + gzip on success.
EXE=/home/baloewe/Overcharged/mpcd/mpcd.out
cd /home/baloewe/Overcharged/runs/Kcalib_lambda
for d in $(cat jobs_susc.txt); do
  ( cd $d; echo "$(date '+%F %T') start $d" >> ../queue.log
    nice -n 19 $EXE -i input.json -o . > log.txt 2>&1 & echo $! > pid; wait $!
    if grep -q "Simulation successful" log.txt; then gzip -f directorfield.dat defects.dat flowfield.dat log.txt 2>/dev/null; touch DONE; fi
    echo "$(date '+%F %T') done  $d" >> ../queue.log ) &
done
wait; echo "$(date '+%F %T') all finished (susc)" >> queue.log
