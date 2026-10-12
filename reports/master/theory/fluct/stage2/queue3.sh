#!/bin/bash
# start the realistic-disorder chunks once the first-round solvers (e12, e15, e20) have finished
cd "$(dirname "$0")"
while pgrep -f "stage2_modes.wl e1[25]|stage2_modes.wl e20" > /dev/null; do sleep 60; done
for c in 00 01 02 03; do setsid -f nohup wolframscript -file run_limited.wl 3_c$c > log3_c$c.txt 2>&1 < /dev/null; sleep 3; done
