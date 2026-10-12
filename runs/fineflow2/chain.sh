#!/bin/bash
# start fineflow2 when fineflow has finished
until [ $(ls /home/baloewe/Physics/Overcharged/runs/fineflow/*/DONE 2>/dev/null | wc -l) -ge 4 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/fineflow2 && ./run_queue.sh jobs.txt 8 > queue.log 2>&1
