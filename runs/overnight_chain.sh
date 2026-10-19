#!/bin/bash
# 2026-10-16: start the overnight sets when the crowded-lambda queue (pid 279457) and the swirl_shear runs have finished.
cd /home/baloewe/Physics/Overcharged/runs
while kill -0 279457 2>/dev/null || ps -eo args | grep -q "^[^ ]*mpcd_swirl.out"; do sleep 300; done
echo "$(date) previous runs done; starting channel_fix (30) and lambda_ellipse (8)"
(cd channel_fix && ./run_queue.sh jobs.txt 14 > queue.log 2>&1) &
(cd lambda_ellipse && ./run_queue.sh jobs.txt 8 > queue.log 2>&1) &
wait; echo "$(date) overnight sets done"
