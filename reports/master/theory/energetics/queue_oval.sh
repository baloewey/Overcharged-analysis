#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
TAG=_oval python3 selfflow.py /home/baloewe/Physics/Overcharged/runs/ovalfine/ba* > selfflow_oval.log 2>&1 &
python3 jetmap.py oval /home/baloewe/Physics/Overcharged/runs/ovalfine/ba* > jetmap_oval.log 2>&1 &
python3 vortdef.py oval /home/baloewe/Physics/Overcharged/runs/ovalfine/ba* > vortdef_oval.log 2>&1 &
SIGMA=2.5 python3 qcrit.py ovals25 /home/baloewe/Physics/Overcharged/runs/ovalfine/ba* > qcrit_oval.log 2>&1 &
wait; echo "$(date) oval analyses done" >> queue_oval.log
