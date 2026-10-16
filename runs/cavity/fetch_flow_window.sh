#!/bin/bash
# Fetch full flow fields of the 10-seed circle window (k 2-10, R 0.8-2.2 l_a; 506 runs, ~2.4 GB) for the defect-centred flow analysis (2026-10-15).
export SSH_AUTH_SOCK=/run/user/1000/gcr/ssh
cd /home/baloewe/Physics/Overcharged/runs/cavity
rsync -a --files-from=flow_window_files.txt baloewe@146.155.46.150:Overcharged/runs/cavity/ server_main/ && echo "$(date) fetched $(wc -l < flow_window_files.txt)"
