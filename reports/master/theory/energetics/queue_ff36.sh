#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
python3 forcefit.py ff36 /home/baloewe/Physics/Overcharged/runs/fineflow/eps* /home/baloewe/Physics/Overcharged/runs/fineflow2/eps* /home/baloewe/Physics/Overcharged/runs/fineflow3/eps* > forcefit_ff36.log 2>&1
python3 vortdef.py ff36 /home/baloewe/Physics/Overcharged/runs/fineflow/eps* /home/baloewe/Physics/Overcharged/runs/fineflow2/eps* /home/baloewe/Physics/Overcharged/runs/fineflow3/eps* > vortdef_ff36.log 2>&1
SIGMA=2.5 python3 qcrit.py ff36s25 /home/baloewe/Physics/Overcharged/runs/fineflow/eps* /home/baloewe/Physics/Overcharged/runs/fineflow2/eps* /home/baloewe/Physics/Overcharged/runs/fineflow3/eps* > qcrit_ff36.log 2>&1
TAG=_ff36 python3 shear_trigger.py /home/baloewe/Physics/Overcharged/runs/fineflow/eps* /home/baloewe/Physics/Overcharged/runs/fineflow2/eps* /home/baloewe/Physics/Overcharged/runs/fineflow3/eps* > shear_trigger_ff36.log 2>&1
echo "$(date) ff36 done" >> queue_ff36.log
