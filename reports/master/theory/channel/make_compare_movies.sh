#!/bin/bash
cd /home/baloewe/Physics/Overcharged/reports/master/theory/channel
for s in sym mixed; do PYTHONPATH=/home/baloewe/Physics/Overcharged/reports/figstyle nice python3 compare_movie.py $s 4000 400 > compare_movie_$s.log 2>&1; done
echo "$(date) movies done"
