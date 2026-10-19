#!/bin/bash
# 2026-10-16: after the overnight sets (overnight_chain.sh pid 281184) start act_robust (12) and channel_curve (24).
cd /home/baloewe/Physics/Overcharged/runs
while kill -0 281184 2>/dev/null; do sleep 300; done
echo "$(date) overnight sets done; starting act_robust and channel_curve"
(cd act_robust && ./run_queue.sh jobs.txt 10 > queue.log 2>&1) &
(cd channel_curve && ./run_queue.sh jobs.txt 12 > queue.log 2>&1) &
wait; echo "$(date) day-2 sets done"
