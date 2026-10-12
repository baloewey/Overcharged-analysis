#!/bin/bash
# Launch every act_* run in directory $1 in parallel (serial code, one core each)
# Usage: ./run_all.sh L100
cd "$(dirname "$0")/$1" || exit 1
EXE=$(realpath ../../../mpcd/mpcd.out)
for d in act_*/; do
    d=${d%/}
    cd $d
    nohup $EXE -i input.json -o . > log.txt 2>&1 &
    echo $! > pid
    cd ..
done
