#!/bin/bash
# switchon2 (2026-10-22): all jobs at once, nice 19. "two NAME TOTAL": passive stage, then switch.py, then active stage. "one NAME TOTAL CKPT": active stage from CKPT.
EXE=/home/baloewe/Overcharged/mpcd/mpcd.out; PY=/home/baloewe/pythonEnv/venv/bin/python; D=/home/baloewe/Overcharged/runs/switchon2
cd $D
run_stage() { ( cd $1; nice -n 19 $EXE -i input.json -o . > log.txt 2>&1 & echo $! > pid; wait $!; grep -q "Simulation successful" log.txt ); }
while read mode name total ck; do
  (
    echo "$(date '+%F %T') start $name" >> queue.log
    if [ "$mode" = two ]; then
      run_stage $name/passive || { echo "$(date '+%F %T') FAILED passive $name" >> queue.log; exit; }
      ck=$D/$name/passive/checkpoint.dat
    fi
    $PY $D/switch.py $ck $D/$name/active/checkpoint_in.dat $total >> queue.log 2>&1
    run_stage $name/active || { echo "$(date '+%F %T') FAILED active $name" >> queue.log; exit; }
    ( cd $name; for f in */directorfield.dat */defects.dat */flowfield.dat */log.txt */avOri.dat; do [ -f $f ] && gzip -f $f; done; touch DONE )
    echo "$(date '+%F %T') done  $name" >> queue.log
  ) &
done < jobs.txt
wait; echo "$(date '+%F %T') all finished" >> queue.log
