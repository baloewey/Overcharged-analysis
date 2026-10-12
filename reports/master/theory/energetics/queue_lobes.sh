#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
until [ $(ls /home/baloewe/Physics/Overcharged/runs/lamact/*/DONE 2>/dev/null | wc -l) -ge 12 ] && grep -q "rotdiff done" queue_rotdiff.log 2>/dev/null; do sleep 120; done
python3 lobes.py 4 > lobes.log 2>&1
echo "$(date) lobes done" >> queue_lobes.log
