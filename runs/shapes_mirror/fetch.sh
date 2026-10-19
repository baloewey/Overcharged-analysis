#!/bin/bash
cd /home/baloewe/Physics/Overcharged/runs/shapes_mirror
rsync -a -e "ssh -o BatchMode=yes" --include='*/' --include='defects.dat.gz' --include='flowfield.dat.gz' --include='input.json' --include='DONE' --include='synopsis.dat' --include='manifest_*.tsv' --exclude='*' baloewe@146.155.46.150:Overcharged/runs/shapes_mirror/ ./ && echo "$(date) fetch done $(find . -name DONE | wc -l)"
