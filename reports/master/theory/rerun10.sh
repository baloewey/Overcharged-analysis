#!/bin/bash
# 10-seed reanalysis chain (27 Oct): defect-based analyses on runs/shapes_server s01-s10 (2250 runs). nice 19, 12 processes.
T=/home/baloewe/Physics/Overcharged/reports/master/theory; export NPROC=12
cd $T/crowd_protect/exchange && nice -n 19 python3 shape_esc2.py > shape_esc2_10.log 2>&1
cd $T/crowd_protect/exchange && nice -n 19 python3 exchange_count2.py 9 > ec2_g9_10.log 2>&1
cd $T/shapes/escape_events && nice -n 19 python3 escape_events2.py > escape_events2_10.log 2>&1
cd $T/crowd_protect/drive_barrier && nice -n 19 python3 drive_barrier.py > drive_barrier_10.log 2>&1 && nice -n 19 python3 drive_test1.py > drive_test1_10.log 2>&1 && nice -n 19 python3 circle_occupancy.py > circle_occupancy_10.log 2>&1
echo "$(date) chain done" >> $T/rerun10.log
