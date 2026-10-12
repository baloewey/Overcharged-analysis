#!/bin/bash
# start the ejection-threshold scan when the negative-winding swirl runs have finished (keeps the machine at <= ~12 MPCD runs)
cd "$(dirname "$0")"
while [ $(ls ../swirl_neg/*/DONE 2>/dev/null | wc -l) -lt 4 ]; do sleep 120; done
B=/home/baloewe/Physics/Overcharged/bin/mpcd_swirl.out
for d in */; do (cd "$d" && setsid -f nohup bash -c "$B -i input.json -o . > log.txt 2>&1; grep -q 'Simulation successful' log.txt && gzip -f flowfield.dat defects.dat; touch DONE" < /dev/null > /dev/null 2>&1); sleep 2; done
echo "$(date '+%F %T') started $(ls -d */ | wc -l) runs"
