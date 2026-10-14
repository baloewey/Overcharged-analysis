#!/bin/bash
# Overnight analyses (2026-10-13): wait for each run set of the after_restart chain, then analyse it (nice 19, one at a time). Log: runs/overnight_analysis.log
R=/home/baloewe/Physics/Overcharged/runs; T=/home/baloewe/Physics/Overcharged/reports/master/theory
waitset() { while [ $(ls $R/$1/*/DONE 2>/dev/null | wc -l) -lt $2 ]; do sleep 600; done; echo "$(date) $1 done"; }
waitset lamsize 18
OUT=tiltlaw_lamsize.txt nice -n 19 python3 $T/energetics/tiltlaw.py $R/lamsize/*/ > $T/energetics/tiltlaw_lamsize.out 2>&1
waitset channel_narrow 16
nice -n 19 python3 $T/channel/channel_an.py $R/channel_narrow/*/ > $T/channel/channel_narrow.txt 2>&1
nice -n 19 python3 $T/channel/channel_current.py $R/channel_narrow/*/ --group > $T/channel/current_narrow.txt 2>&1
waitset channel_qscan 6
nice -n 19 python3 $T/channel/channel_an.py $R/channel_qscan/*/ > $T/channel/channel_qscan.txt 2>&1
nice -n 19 python3 $T/channel/channel_current.py $R/channel_qscan/*/ --group > $T/channel/current_qscan.txt 2>&1
waitset channel_mixed 10
nice -n 19 python3 $T/channel/channel_current.py $R/channel_mixed/*/ --group > $T/channel/current_mixed.txt 2>&1
waitset channel_mixedW 12
nice -n 19 python3 $T/channel/channel_current.py $R/channel_mixedW/*/ --group > $T/channel/current_mixedW.txt 2>&1
echo "$(date) all analyses done"
