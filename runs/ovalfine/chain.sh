#!/bin/bash
# start the oval fine-flow runs once the hysteresis runs are done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/hysteresis/*/DONE 2>/dev/null | wc -l) -ge 16 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/ovalfine && ./run_queue.sh jobs.txt 12 > queue.log 2>&1
