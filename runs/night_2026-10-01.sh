#!/bin/bash
# Overnight analysis, 2026-10-01. Waits for the office runs, then runs the standard analyses.
# Read-only on the runs; writes figures/tables under runs/night_2026-10-01/ and a log next to this script.
# Launch: setsid -f nohup ./night_2026-10-01.sh > night_2026-10-01.log 2>&1 < /dev/null
set -u
cd "$(dirname "$0")"
R=$PWD; O=$R/night_2026-10-01; mkdir -p $O
source $R/../.venv/bin/activate 2>/dev/null
log() { echo "$(date '+%F %T') $*"; }

alldone() {  # $1 = jobs file
    while read d; do [ -f "$d/DONE" ] || return 1; done < "$1"; return 0; }

log "waiting for: cassini P1, P2, P3 and the zeta = 0.09 seeds"
until alldone cassini/jobs_pilot1.txt && alldone cassini/jobs_pilot2.txt && alldone cassini/jobs_pilot3.txt \
      && alldone cavity/jobs_h3_act0.09_seeds.txt; do sleep 120; done
log "all office runs done"

# 1. accumulated-charge figure (picks up the half-length-34 ovals automatically)
( cd cassini && python3 fig_accumulated_charge.py ) && cp cassini/fig_accumulated_charge.png cassini/fig_accumulated_charge.pdf $O/
log "1 accumulated charge figure done"

# 2. small ovals + high k: bulk/edge numbers, lobe profiles, cut-free criteria, label tracking
SMALL=$(ls -d cassini/P3/H34_ba*_k06/s0? | tr '\n' ' ')
HIGHK=$(ls -d cassini/P3/H48_ba1.02_k?? | tr '\n' ' ')
( cd cassini && python3 cassini_charge.py ${SMALL//cassini\//} ${HIGHK//cassini\//} ) > $O/cassini_charge_P3.txt 2>&1
( cd cassini && python3 lobe_profiles.py ${SMALL//cassini\//} ${HIGHK//cassini\//} --out $O/lobe_profiles_P3 ) > $O/lobe_profiles_P3.txt 2>&1
( cd cavity && python3 edge_criteria.py ${SMALL//cassini/..\/cassini} ${HIGHK//cassini/..\/cassini} --out $O/edge_criteria_P3 ) > $O/edge_criteria_P3.txt 2>&1
# label tracking: small ovals grouped by shape, with the circles of comparable size for reference
LABS=""; RUNS=""
for d in cassini/P3/H34_ba*_k06; do ba=${d##*_ba}; ba=${ba%_k06}; for s in $d/s0?; do RUNS="$RUNS ../$s"; LABS="$LABS oval_b/a_${ba}_H1.53la"; done; done
for s in cavity/server_k06/rho0.06/k06/s0? cavity/server_k06/rho0.05/k06/s0?; do r=${s#cavity/}; RUNS="$RUNS $r"; case $s in *rho0.06*) LABS="$LABS circle_R1.43la";; *) LABS="$LABS circle_R1.72la";; esac; done
for d in $HIGHK; do RUNS="$RUNS ../$d"; LABS="$LABS oval_1.02_H2.16la_k${d##*_k}"; done
( cd cavity && python3 bc_labels.py $RUNS --labels $LABS --out $O/bc_labels_P3 ) > $O/bc_labels_P3.txt 2>&1
log "2 small ovals / high k analyses done"

# 3. renders (snapshots + movies) of all finished active ovals not yet rendered
( cd cassini && for d in P1/H* P2/H* P3/H48_ba1.02_k??; do n=$(basename $d); [ -f movies/${n}_last.png ] || python3 render_cassini.py $d --outdir movies; done; \
  for d in P3/H34_ba*_k06; do python3 render_cassini.py $d/s01 --outdir movies/P3_s01; done ) > $O/render.txt 2>&1
( cd cassini && python3 lobe_rotation.py P1/H* P2/H* P3/H34_ba*_k06/s0? P3/H48_ba1.02_k?? ) > $O/lobe_rotation_all.txt 2>&1
log "3 renders and lobe rotation done"

# 4. zeta = 0.09 with 3 seeds: site rule, circulation, enclosed charge
( cd cavity && R09="H3pilot/act0.09/k0[2-9] H3pilot/act0.09/k10 H3pilot_seeds/act0.09/s0[23]/k0[2-9] H3pilot_seeds/act0.09/s0[23]/k10" && \
  python3 edge_sites.py $R09 --no-gr --out $O/act009_3seeds_sites > $O/act009_3seeds_sites.txt 2>&1; \
  python3 circulation.py $R09 --out $O/circ_act009_3seeds.csv > $O/circ_act009_3seeds.txt 2>&1 )
log "4 zeta = 0.09 analyses done"

# 5. server progress (read-only)
timeout 30 ssh -o BatchMode=yes baloewe@146.155.46.150 'cd ~/Overcharged/runs/cavity; date; echo "batch 2+3 done $(grep -c " done " queue_batch23.log)/1304 fail $(grep -c FAIL queue_batch23.log)"; pgrep -c mpcd; df -h / | tail -1' > $O/server_status.txt 2>&1
log "5 server status written; finished"
