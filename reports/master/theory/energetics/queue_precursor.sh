#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
while pgrep -f "steadyE.py oval" > /dev/null; do sleep 60; done
python3 precursor.py 3 > precursor.log 2>&1
echo "$(date) precursor done" >> queue_precursor.log
