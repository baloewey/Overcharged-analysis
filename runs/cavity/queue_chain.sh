#!/bin/bash
# Chain batches on the server: once batch 1's queue has STARTED all its jobs (its long runs may
# still be going), launch one pool of 96 running batch 2 then batch 3 (jobs_batch23.txt).
# Brief oversubscription while batch 1's tail finishes is intentional (keeps all threads busy).
cd "$(dirname "$0")"
started() { grep -cE " start | skip " "$1" 2>/dev/null || echo 0; }
N1=$(wc -l < jobs_batch1.txt)
cat jobs_batch2.txt jobs_batch3.txt > jobs_batch23.txt
echo "$(date '+%F %T') chain: waiting for batch 1 to start all $N1 jobs"
while [ "$(started queue_batch1.log)" -lt "$N1" ]; do sleep 120; done
echo "$(date '+%F %T') batch 1 fully started ($(pgrep -c mpcd.out) solvers running); launching batches 2+3 ($(wc -l < jobs_batch23.txt) jobs, 96 parallel)"
./run_queue.sh jobs_batch23.txt 96 > queue_batch23.log 2>&1
echo "$(date '+%F %T') chain: batches 2+3 queue finished"
