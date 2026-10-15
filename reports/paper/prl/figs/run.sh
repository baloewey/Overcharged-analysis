#!/bin/bash
# Run a PRL figure script with the house capitalisation but no font scaling (figures are drawn at printed size).
cd "$(dirname "$0")" && FIGSTYLE_SCALE=1 PYTHONPATH=/home/baloewe/Physics/Overcharged/reports/figstyle:$(pwd)${PYTHONPATH:+:$PYTHONPATH} exec python3 "$@"
