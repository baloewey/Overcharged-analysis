#!/bin/bash
# Wait for the L100 sweep to finish, then launch the L200 finite-size checks.
cd "$(dirname "$0")"
PIDS=$(cat L100/act_*/pid)
echo "$(date '+%F %T') waiting for L100 pids: $PIDS"
while :; do
    alive=0
    for p in $PIDS; do kill -0 $p 2>/dev/null && alive=1; done
    [ $alive -eq 0 ] && break
    sleep 60
done
echo "$(date '+%F %T') L100 finished; launching L200"
./run_all.sh L200
echo "$(date '+%F %T') launched: $(cat L200/act_*/pid | tr '\n' ' ')"
