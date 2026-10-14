#!/bin/bash
# Run a figure script with the house style (capitalised labels, larger fonts): reports/figstyle/run.sh script.py [args]
PYTHONPATH=/home/baloewe/Physics/Overcharged/reports/figstyle${PYTHONPATH:+:$PYTHONPATH} exec python3 "$@"
