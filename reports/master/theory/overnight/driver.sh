#!/bin/bash
# Overnight driver (detached): analyses each batch as soon as it is complete, writes summaries next to this script.
cd "$(dirname "$0")"
R=/home/baloewe/Physics/Overcharged/runs
wait_batch() { while true; do n=$(ls -d $R/$1/*/ 2>/dev/null | wc -l); d=$(ls $R/$1/*/DONE 2>/dev/null | wc -l); [ "$n" -gt 0 ] && [ "$d" -ge "$n" ] && break; sleep 120; done; }
echo "$(date '+%F %T') driver start"
wait_batch swirl_neg;       python3 counts.py swirl_neg       > swirl_neg_summary.txt 2>/dev/null;       echo "$(date '+%F %T') swirl_neg analysed"
wait_batch swirl_ellipse;   python3 counts.py swirl_ellipse   > swirl_ellipse_summary.txt 2>/dev/null;   echo "$(date '+%F %T') swirl_ellipse analysed"
while [ ! -s $R/swirl_threshold/queue.log ]; do sleep 120; done
wait_batch swirl_threshold; python3 counts.py swirl_threshold > swirl_threshold_summary.txt 2>/dev/null; echo "$(date '+%F %T') swirl_threshold analysed"
echo "$(date '+%F %T') driver done"
