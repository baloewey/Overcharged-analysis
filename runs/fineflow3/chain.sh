#!/bin/bash
# start fineflow3 once the bulk-only runs are done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/bulkfine/*/DONE 2>/dev/null | wc -l) -ge 10 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/runs/fineflow3 && ./run_queue.sh jobs.txt 20 > queue.log 2>&1
