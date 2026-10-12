#!/bin/bash
# tonight (user: max 8 cores / 16 threads): neglayer -> lowlam_fine (12 parallel) -> lamact (12 parallel)
until [ $(ls /home/baloewe/Physics/Overcharged/runs/neglayer/*/DONE 2>/dev/null | wc -l) -ge 14 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/lowlam_fine && ./run_queue.sh jobs.txt 12 > queue.log 2>&1
cd /home/baloewe/Physics/Overcharged/runs/lamact && ./run_queue.sh jobs.txt 12 > queue.log 2>&1
