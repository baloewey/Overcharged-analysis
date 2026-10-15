#!/bin/bash
# After the 10-seed sweep tables (2026-10-15): flow/tilt table now, then (once rerun_seeds10.sh has exited) the plots and report figures that read the tables.
cd /home/baloewe/Physics/Overcharged/runs/cavity
echo "$(date) flow_tilt start"; nice -n 19 python3 flow_tilt.py --workers 6 > logs_seeds10_flowtilt.log 2>&1; echo "$(date) flow_tilt $?"
P=$(pgrep -f rerun_seeds10.sh | head -1); while [ -n "$P" ] && kill -0 $P 2>/dev/null; do sleep 60; done
echo "$(date) tables done"
FS=/home/baloewe/Physics/Overcharged/reports/figstyle/run.sh
for s in plots_criteria.py plots_flow.py plots_fourier.py; do nice -n 19 $FS $s > logs_seeds10_$s.log 2>&1; echo "$(date) $s $?"; done
M=/home/baloewe/Physics/Overcharged/reports/master
for s in $M/scripts/fig_onset_phase.py $M/theory/effbc/spacing_criterion.py $M/theory/stability/tilt_onset_k.py; do (cd $(dirname $s) && nice -n 19 $FS $(basename $s) > /home/baloewe/Physics/Overcharged/runs/cavity/logs_seeds10_$(basename $s).log 2>&1); echo "$(date) $(basename $s) $?"; done
echo "$(date) post done"
