#!/bin/bash
# Start the crowding runs after the channel_current queue (pid 115171) finishes.
cd /home/baloewe/Physics/Overcharged/runs/channel_crowd
while kill -0 115171 2>/dev/null; do sleep 300; done
echo "$(date) channel_current finished; starting channel_crowd"
./run_queue.sh jobs.txt 22 > queue.log 2>&1
echo "$(date) channel_crowd done"
