#!/bin/bash
# Wait for the L100 active-length sweep to finish, then launch the fixed-R cavity scan.
cd "$(dirname "$0")"
PIDS=$(cat ../activeLength/L100/act_*/pid)
echo "$(date '+%F %T') waiting for L100 pids: $(echo $PIDS)"
while :; do
    alive=0
    for p in $PIDS; do kill -0 $p 2>/dev/null && alive=1; done
    [ $alive -eq 0 ] && break
    sleep 60
done
EXE=$(realpath ../../mpcd/mpcd.out)
for d in Rscan_R31.83/k*/; do
    cd $d
    nohup $EXE -i input.json -o . > log.txt 2>&1 &
    echo $! > pid
    cd - > /dev/null
done
echo "$(date '+%F %T') launched: $(cat Rscan_R31.83/k*/pid | tr '\n' ' ')"
