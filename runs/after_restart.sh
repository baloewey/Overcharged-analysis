#!/bin/bash
# Start the runs held for the office-PC restart (2026-10-13). Thread budget 16:
#   2 long channels (~20 h)  +  14 slots: lambda-vs-size circles (18 runs, ~3-4 h), then narrow channels (16 runs, ~2 h).
# Usage (after the reboot): ~/Physics/Overcharged/runs/after_restart.sh
cd /home/baloewe/Physics/Overcharged/runs/channel_long && setsid nohup ./run_queue.sh jobs.txt 2 > queue.log 2>&1 < /dev/null &
cd /home/baloewe/Physics/Overcharged/runs && setsid nohup bash -c 'cd lamsize && ./run_queue.sh jobs.txt 14 > queue.log 2>&1; cd ../channel_narrow && ./run_queue.sh jobs.txt 14 > queue.log 2>&1' > /dev/null 2>&1 < /dev/null &
echo "started channel_long (2) and lamsize (18) -> channel_narrow (16), 14 parallel"
