#!/bin/bash
# Chain (2026-10-15): start the shape 10-seed extension (+k then -k) when the mirror shape sweep queue has finished.
LOG=~/Overcharged/runs/chain_shapes_seeds10.log
Q=1087080   # run_queue.sh jobs_shape_mirror.txt (checked 2026-10-15)
echo "$(date) waiting for the mirror shape queue (pid $Q)" >> $LOG
while kill -0 $Q 2>/dev/null; do sleep 300; done
echo "$(date) mirror queue finished: $(grep -c "done  " ~/Overcharged/runs/shapes_mirror/queue_shape_mirror.log) done lines; starting shape seeds" >> $LOG
cd ~/Overcharged/runs && ./shapes_mirror/run_queue.sh jobs_shapes_seeds10.txt 96 > queue_shapes_seeds10.log 2>&1
echo "$(date) shape seeds queue exited" >> $LOG
