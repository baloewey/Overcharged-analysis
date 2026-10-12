#!/bin/bash
# H3 test inside the turbulent regime: activities 0.06 and 0.11 at R = 1.434 l_a (same R/l_a as the
# pilots at 0.05/0.07/0.09).
#   Stage 1: l_a in a 200^2 periodic box (same inputs as runs/activeLength/L200, 2 seeds per activity)
#   Stage 2: l_a from active_length.py (seed mean), cavities k = 2..10 x 3 seeds per activity,
#            inputs identical to H3pilot_seeds apart from act, R, seed; legacy tracker
#   Stage 3: run the 54 cavities, 18 in parallel (run_queue.sh)
# Restart-safe: finished stage-1 runs (log with "Simulation successful") are not rerun.
# Launch: setsid -f nohup ./chain_h3_high.sh > chain_h3_high.log 2>&1 < /dev/null
set -u
cd "$(dirname "$0")"
HERE=$PWD
AL=$(realpath ../activeLength)
EXE=$(realpath ../../mpcd/mpcd.out)
ACTS="0.060 0.110"
log() { echo "$(date '+%F %T') $*"; }

# ---- stage 1: l_a runs ----
pids=""
for a in $ACTS; do
    for s in 1 2; do
        d=$AL/L200/act_${a}_s$s
        mkdir -p $d
        if [ -f $d/log.txt ] && grep -q "Simulation successful" $d/log.txt; then log "skip (done) $d"; continue; fi
        seed=$(python3 -c "print(20000 + int(round(float('$a')*1000))*10 + $s)")
        cat > $d/input.json <<EOF
{
    "domain":           [200, 200],
    "simSteps":         25000,
    "warmUp":           0,
    "seed":             $seed,
    "collOp":           21,
    "lc":               1,
    "species":          [ { "dens": 20, "act": $a } ],
    "domainWalls":      1,
    "defectsOut":       20,
    "defectTracker":    0,
    "checkpointOut":    5000,
    "synopsisOut":      1
}
EOF
        ( cd $d && exec $EXE -i input.json -o . > log.txt 2>&1 ) &
        echo $! > $d/pid
        pids="$pids $!"
        log "start l_a run $d (seed $seed)"
    done
done
for p in $pids; do wait $p; done
for a in $ACTS; do for s in 1 2; do
    d=$AL/L200/act_${a}_s$s
    grep -q "Simulation successful" $d/log.txt || { log "FAIL $d, stopping"; exit 1; }
done; done
log "stage 1 done"

# ---- stage 2: l_a and cavity inputs ----
cd $AL
python3 active_length.py L200/act_0.060_s1 L200/act_0.060_s2 L200/act_0.110_s1 L200/act_0.110_s2 \
    --L 200 200 --out L200/active_length_h3high.csv > L200/active_length_h3high.txt 2>&1 || { log "active_length.py failed"; exit 1; }
cd $HERE
python3 - <<'EOF' || { echo "cavity generation failed"; exit 1; }
import csv, json, os, subprocess, collections, numpy as np
rows = [r for r in csv.DictReader(open("../activeLength/L200/active_length_h3high.csv")) if r["estimator"] == "all"]
by = collections.defaultdict(list)
for r in rows:
    by[round(float(r["act"]), 3)].append((float(r["l_active"]), float(r["l_active_err"]), r["status"]))
jobs = open("jobs_h3_high.txt", "w")
import os; os.makedirs("H3high", exist_ok=True)
with open("H3high/l_a_used.txt", "w") as fo:
    fo.write("act  l_a(seed1)  l_a(seed2)  l_a(mean)  err  R = 1.434 l_a\n")
    for act, v in sorted(by.items()):
        la = np.array([x[0] for x in v]); er = np.array([x[1] for x in v])
        w = 1 / er**2; m = (la * w).sum() / w.sum(); e = 1 / np.sqrt(w.sum())
        R = round(1.434 * m, 4)
        fo.write(f"{act:.3f}  " + "  ".join(f"{x:.2f}±{y:.2f} ({st})" for x, y, st in v) + f"  {m:.2f}  {e:.2f}  {R}\n")
        for s in (1, 2, 3):
            for k in range(2, 11):
                d = f"H3high/act{act:.2f}/s{s:02d}/k{k:02d}"
                seed = int(round(act * 1000)) * 1000 + 100 * s + k
                subprocess.run(["python3", "gen_cavity.py", d, "--k", str(k), "--R", str(R), "--act", str(act),
                                "--seed", str(seed), "--warmup", "10000", "--steps", "40000", "--defout", "20",
                                "--dirout", "500", "--flowout", "500", "--checkpoint", "10000"],
                               check=True, stdout=subprocess.DEVNULL)
                p = f"{d}/input.json"; j = json.load(open(p)); j["defectTracker"] = 0
                json.dump(j, open(p, "w"), indent=4)
                jobs.write(os.path.abspath(d) + "\n")
jobs.close()
EOF
log "stage 2 done: $(wc -l < jobs_h3_high.txt) cavity runs; l_a:"
cat H3high/l_a_used.txt

# ---- stage 3: cavities ----
./run_queue.sh jobs_h3_high.txt 18 > queue_h3_high.log 2>&1
log "stage 3 done: $(grep -c ' done ' queue_h3_high.log) done, $(grep -c FAIL queue_h3_high.log) failed"
