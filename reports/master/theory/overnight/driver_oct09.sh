#!/bin/bash
# Detached fallback (2026-10-09 evening): when the runs finish, write basic summaries even if the Claude session is gone.
R=/home/baloewe/Physics/Overcharged/runs
until [ $(ls $R/tiltdiag/*/DONE 2>/dev/null | wc -l) -ge 6 ] && [ $(ls $R/rotvisc/*/DONE 2>/dev/null | wc -l) -ge 6 ]; do sleep 120; done
cd /home/baloewe/Physics/Overcharged/reports/master/theory/overnight && python3 counts.py tiltdiag > tiltdiag_summary.txt 2>&1
cd $R/Kcalib && python3 rot_visc.py $R/rotvisc/modes_L64_s1 $R/rotvisc/modes_L64_s2 $R/rotvisc/modes_L64_s3 $R/rotvisc/modes_L64_s4 > $R/rotvisc/rot_visc_modes.txt 2>&1
echo "$(date) driver done" >> /home/baloewe/Physics/Overcharged/reports/master/theory/overnight/driver_oct09.log
