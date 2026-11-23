#!/bin/bash
# Run every job in jobs.txt in parallel (nice 19) with the cassini binary (as the shape sweep); gzip outputs; DONE marker.
D=$(cd "$(dirname "$0")" && pwd); EXE=$(realpath ~/Overcharged/niamh-overcharged-cassini/mpcd.out); cd $D
while read n; do
  ( cd $n; [ -f DONE ] && exit 0; echo "$(date '+%F %T') start $n" >> ../queue.log
    nice -n 19 $EXE -i input.json -o . > log.txt 2>&1 & echo $! > pid; wait $!
    if grep -q "Simulation successful" log.txt; then gzip -f defects.dat directorfield.dat flowfield.dat log.txt 2>/dev/null; rm -f checkpoint.dat; date > DONE; echo "$(date '+%F %T') done  $n" >> ../queue.log
    else echo "$(date '+%F %T') FAIL  $n" >> ../queue.log; fi ) &
done < ${JOBS:-jobs.txt}
wait; echo "$(date '+%F %T') all finished" >> queue.log
