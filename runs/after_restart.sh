#!/bin/bash
# Start the runs held for the office-PC restart (2026-10-13). Thread budget 16:
#   2 long channels (~20 h)  +  14 slots: lambda-vs-size circles (18 runs, ~3-4 h), then narrow channels (16 runs, ~2 h), then the channel q-scan (6 runs, ~2-3 h), then negative/mixed/unequal-winding channels (10 runs, ~3 h), then mixed channels at W 10/16/24/60 (12 runs, ~3 h).
# Usage (after the reboot): ~/Physics/Overcharged/runs/after_restart.sh
cd /home/baloewe/Physics/Overcharged/runs/channel_long && setsid nohup ./run_queue.sh jobs.txt 2 > queue.log 2>&1 < /dev/null &
cd /home/baloewe/Physics/Overcharged/runs && setsid nohup bash -c 'cd lamsize && ./run_queue.sh jobs.txt 14 > queue.log 2>&1; cd ../channel_narrow && ./run_queue.sh jobs.txt 14 > queue.log 2>&1; cd ../channel_qscan && ./run_queue.sh jobs.txt 14 > queue.log 2>&1; cd ../channel_mixed && ./run_queue.sh jobs.txt 14 > queue.log 2>&1; cd ../channel_mixedW && ./run_queue.sh jobs.txt 14 > queue.log 2>&1' > /dev/null 2>&1 < /dev/null &
echo "started channel_long (2) and lamsize (18) -> channel_narrow (16) -> channel_qscan (6) -> channel_mixed (10) -> channel_mixedW (12), 14 parallel"
