#!/bin/bash
cd /home/baloewe/Physics/Overcharged/runs/negseed2
until grep -q "queue finished" queue.log; do sleep 120; done
/home/baloewe/Physics/Overcharged/reports/figstyle/run.sh analyse.py > analyse.log 2>&1; echo "$(date) analysed" >> analyse.log
