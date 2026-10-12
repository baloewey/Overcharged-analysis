#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
until [ $(ls /home/baloewe/Physics/Overcharged/runs/bulkfine/*/DONE 2>/dev/null | wc -l) -ge 10 ]; do sleep 120; done
python3 forcefit.py bulk /home/baloewe/Physics/Overcharged/runs/bulkfine/eps* > forcefit_bulk.log 2>&1
echo "$(date) forcefit bulk done" >> queue_bulk.log
