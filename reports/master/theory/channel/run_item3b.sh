#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel
until grep -q "batch1 done" run_item3.log; do sleep 30; done
nice python3 capacitor.py > capacitor.log 2>&1
nice python3 pressure.py > pressure.log 2>&1
echo "$(date) item3 batch2 done"
