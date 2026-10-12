#!/bin/bash
# local escape analyses on the bulk-only fine-output runs once they finish (selfflow, shear_trigger, energetics_ts)
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
until [ $(ls /home/baloewe/Physics/Overcharged/runs/bulkfine/*/DONE 2>/dev/null | wc -l) -ge 10 ]; do sleep 120; done
R=/home/baloewe/Physics/Overcharged/runs/bulkfine
TAG=_bulk python3 selfflow.py $R/eps* > selfflow_bulk.log 2>&1 &
TAG=_bulk python3 shear_trigger.py $R/eps* > shear_trigger_bulk.log 2>&1 &
python3 energetics_ts.py $R/eps* bulk > ets_bulk.log 2>&1 &
wait; echo "$(date) done" >> queue_bulk.log
