#!/bin/bash
# Fetch defects + inputs of the 10-seed circle production (s04-s10, + and - winding) into server_main (2026-10-15).
export SSH_AUTH_SOCK=/run/user/1000/gcr/ssh
cd /home/baloewe/Physics/Overcharged/runs/cavity/server_main
for set in main_act0.07 main_act0.07_neg; do
  rsync -a --prune-empty-dirs --include='*/' --include='s0[4-9]/input.json' --include='s0[4-9]/defects.dat.gz' --include='s0[4-9]/DONE' \
        --include='s10/input.json' --include='s10/defects.dat.gz' --include='s10/DONE' --exclude='*' \
        baloewe@146.155.46.150:Overcharged/runs/cavity/$set/ ./$set/
  echo "$(date) $set fetched"
done
