#!/bin/bash
# After the flow fields of the 10-seed circle window are fetched: defect-centred flow and wall-band shear on 10 seeds (2026-10-15).
cd /home/baloewe/Physics/Overcharged/reports/master/theory/tiltcrowd
until grep -q fetched /home/baloewe/Physics/Overcharged/runs/cavity/fetch_flow_window.log 2>/dev/null; do sleep 60; done
FS=/home/baloewe/Physics/Overcharged/reports/figstyle/run.sh
MAIN10=1 TAG=_main10 nice -n 19 $FS defect_flow.py 6 > defect_flow_main10.log 2>&1; echo "$(date) defect_flow $?"
MAIN10=1 TAG=_main10 nice -n 19 $FS wallshear.py 6 > wallshear_main10.log 2>&1; echo "$(date) wallshear $?"
