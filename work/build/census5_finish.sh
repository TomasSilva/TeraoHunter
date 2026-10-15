#!/bin/sh
# post-process the five-pencil census: classes, exact lifting, explicit solutions
cd /Users/agent/Projects/TeraoHunter/work/build
PY=/Users/agent/Projects/TeraoHunter/.venv/bin/python
$PY census5_post.py > census5_classes.txt 2>&1
cd ../tc15
NP=4
for part in $(seq 0 $((NP-1))); do $PY exact_tc_lattice.py 17 8 ../build/census5_classes.jsonl $part $NP > ../build/census5_lift_$part.log 2>&1 & done
wait
cd ../build
cat census5_classes_tc_*.out | grep -v "^$" > census5_lift.txt
echo FINISHED
