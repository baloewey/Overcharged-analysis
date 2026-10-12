#!/bin/bash
R=/home/baloewe/Physics/Overcharged/runs
until [ $(ls $R/actlayout/*/DONE 2>/dev/null | wc -l) -ge 8 ]; do sleep 120; done
python3 energetics.py $R/actlayout/eps* > energetics_actlayout.log 2>&1
until [ $(ls $R/fineflow/*/DONE 2>/dev/null | wc -l) -ge 4 ]; do sleep 120; done
python3 energetics.py $R/fineflow/eps* > energetics_fineflow.log 2>&1
echo "$(date) done" >> queue_after.log
