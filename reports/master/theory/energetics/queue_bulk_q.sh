#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
until [ $(ls /home/baloewe/Physics/Overcharged/runs/bulkfine/*/DONE 2>/dev/null | wc -l) -ge 10 ]; do sleep 120; done
SIGMA=2.5 python3 qcrit.py bulks25 /home/baloewe/Physics/Overcharged/runs/bulkfine/eps* > qcrit_bulk.log 2>&1
echo "$(date) qcrit bulk done" >> queue_bulk.log
