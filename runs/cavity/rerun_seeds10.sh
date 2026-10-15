#!/bin/bash
# Rerun the circle pipeline on the 10-seed sweep (2026-10-15). Defect-based steps; flow steps follow once flow_red is fetched.
cd /home/baloewe/Physics/Overcharged/runs/cavity
P=python3
echo "$(date) start"
nice -n 19 $P criteria_all.py --workers 10 > logs_seeds10_criteria.log 2>&1; echo "$(date) criteria_all $?"
nice -n 19 $P wall_fourier_all.py --workers 10 > logs_seeds10_fourier.log 2>&1; echo "$(date) wall_fourier_all $?"
nice -n 19 $P mirror_analysis.py > logs_seeds10_mirror.log 2>&1; echo "$(date) mirror $?"
nice -n 19 $P main_sweep_analysis.py > logs_seeds10_main.log 2>&1; echo "$(date) main_sweep $?"
