#!/bin/bash
# start the seeded hysteresis runs once the circle lambda scan is done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/lambda_circle/*/DONE 2>/dev/null | wc -l) -ge 33 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/hysteresis && ./run_queue.sh jobs.txt 16 > queue.log 2>&1
