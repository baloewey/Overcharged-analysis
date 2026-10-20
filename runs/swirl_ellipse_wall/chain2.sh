#!/bin/bash
cd /home/baloewe/Physics/Overcharged/runs/swirl_ellipse_wall
until [ $(for d in $(cat jobs.txt); do [ -f $d/DONE ] && echo x; done | wc -l) -ge 8 ]; do sleep 60; done
./run_queue.sh jobs2.txt 8 > queue2.log 2>&1
