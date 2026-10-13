#!/bin/bash
# after tonight's lamact runs (user limit 16 threads: 8 channel runs + lobes analysis 4)
until [ $(ls /home/baloewe/Physics/Overcharged/runs/lamact/*/DONE 2>/dev/null | wc -l) -ge 12 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/channel && ./run_queue.sh jobs.txt 8 > queue.log 2>&1
