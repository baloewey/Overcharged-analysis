#!/bin/bash
# two long channels alongside the L192 scan once lamact is done (user limit 16 threads: 8 + 2 + lobes 4)
until [ $(ls /home/baloewe/Physics/Overcharged/runs/lamact/*/DONE 2>/dev/null | wc -l) -ge 12 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/channel_long && ./run_queue.sh jobs.txt 2 > queue.log 2>&1
