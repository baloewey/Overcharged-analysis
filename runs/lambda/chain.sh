#!/bin/bash
# start the lambda runs once fineflow3 is done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/fineflow3/*/DONE 2>/dev/null | wc -l) -ge 24 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/lambda && ./run_queue.sh jobs.txt 20 > queue.log 2>&1
