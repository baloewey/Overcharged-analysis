#!/bin/bash
while kill -0 536431 2>/dev/null; do sleep 20; done
cd /home/baloewe/Physics/Overcharged/reports/master/theory/shapes
echo "$(date) fetch: $(find ../../../../runs/shapes_mirror -name DONE | wc -l) DONE"
nice -n 5 python3 mirror_shape_analysis.py && echo "$(date) mirror analysis done"
