#!/bin/bash
# start the circle lambda scan (outline plan A) once runs/lambda is done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/lambda/*/DONE 2>/dev/null | wc -l) -ge 20 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/lambda_circle && ./run_queue.sh jobs.txt 28 > queue.log 2>&1
