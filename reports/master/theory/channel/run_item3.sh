#!/bin/bash
# analysis item 3 (channels), 2026-10-19, detached
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel
RDIR=/home/baloewe/Physics/Overcharged/runs/channel_fix PAT="crowd_*" OUT=crowd_an_fix.txt nice python3 crowd_an.py > crowd_an_fix.log 2>&1 &
nice python3 current_stats.py > current_stats.log 2>&1 &
nice python3 net_flow.py > net_flow.log 2>&1 &
wait; echo "$(date) item3 batch1 done"
