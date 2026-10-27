#!/bin/bash
# Passive imposed-vortex controls, server: all jobs at once, nice 19, swirl binary (anchmod branch + body force, bin/mpcd_swirl.out).
cd "$(dirname "$0")"; B=$HOME/Overcharged/bin/mpcd_swirl.out
for d in $(cat jobs.txt); do (cd "$d" && setsid -f nohup nice -n 19 bash -c "$B -i input.json -o . > log.txt 2>&1; grep -q 'Simulation successful' log.txt && gzip -f flowfield.dat defects.dat directorfield.dat; date > DONE" < /dev/null > /dev/null 2>&1); sleep 1; done
echo "$(date '+%F %T') started $(wc -l < jobs.txt) runs"
