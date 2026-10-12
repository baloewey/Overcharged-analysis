#!/bin/bash
# Run every directory listed in a jobs file, NPAR at a time (serial code, one core each).
# Jobs are started in file order (gen_sweep.py writes most expensive first).
#
# When a run finishes successfully:
#   - a DONE marker is written,
#   - large outputs are gzip-compressed (directorfield, flowfield, defects, topochargefield,
#     checkpoint, log.txt) to save SSD space. The analysis scripts read .gz transparently.
# Runs with a DONE marker (or an uncompressed log with "Simulation successful") are skipped,
# so the queue can be stopped and restarted safely. Runs killed mid-way restart from scratch.
#
# Usage:   ./run_queue.sh JOBS.txt [NPAR=12]
# Detach:  setsid nohup ./run_queue.sh JOBS.txt 96 > queue.log 2>&1 < /dev/null &
# Stop:    pkill -f run_queue.sh; then kill the running solvers: kill $(cat .../s*/pid)

JOBS=${1:?usage: run_queue.sh jobs.txt [NPAR]}
NPAR=${2:-12}
EXE=$(realpath /home/baloewe/Physics/Overcharged/bin/mpcd_wallact2.out)
export EXE

run_one() {
    d=$1
    if [ -f "$d/DONE" ] || { [ -f "$d/log.txt" ] && grep -q "Simulation successful" "$d/log.txt"; }; then
        echo "$(date '+%F %T') skip (done) $d"
        return
    fi
    cd "$d" || return
    echo "$(date '+%F %T') start $d"
    "$EXE" -i input.json -o . > log.txt 2>&1 &
    echo $! > pid
    wait $!
    if grep -q "Simulation successful" log.txt; then
        echo "$(date '+%F %T') done $d" > DONE
        for f in directorfield.dat flowfield.dat defects.dat topochargefield.dat checkpoint.dat log.txt; do
            [ -f "$f" ] && gzip -6 -f "$f"
        done
        echo "$(date '+%F %T') done  $d"
    else
        echo "$(date '+%F %T') FAIL  $d (see log.txt)"
    fi
}
export -f run_one

echo "$(date '+%F %T') queue start: $(wc -l < "$JOBS") jobs, $NPAR in parallel"
xargs -a "$JOBS" -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
echo "$(date '+%F %T') queue finished"
