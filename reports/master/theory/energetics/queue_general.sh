#!/bin/bash
# shape-sweep energetics (ellipses + circles) once the activity-layout runs are done
until [ $(ls /home/baloewe/Physics/Overcharged/runs/actlayout/*/DONE 2>/dev/null | wc -l) -ge 8 ]; do sleep 120; done
SW=/home/baloewe/Physics/Overcharged/runs/shapes_server
python3 energetics_general.py shapes $(ls -d $SW/ellipse/*/*/*/s0* $SW/circle/circle/*/*/s0*) > etsg_shapes.log 2>&1
echo "$(date) done" >> queue_general.log
