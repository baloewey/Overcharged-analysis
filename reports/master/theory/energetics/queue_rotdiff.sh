#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
until grep -q "precursor done" queue_precursor.log 2>/dev/null; do sleep 60; done
python3 rotdiff.py 3 > rotdiff.log 2>&1
echo "$(date) rotdiff done" >> queue_rotdiff.log
