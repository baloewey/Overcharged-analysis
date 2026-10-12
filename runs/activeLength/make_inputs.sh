#!/bin/bash
# Active length-scale sweep: LxL PBC bulk active nematic, act = 0.050 ... 0.090
# Only defect data (plus a checkpoint for extending runs) is written.
# Usage: ./make_inputs.sh L STEPS ["act1 act2 ..."]   -> creates L${L}/act_*/input.json
L=${1:-100}
STEPS=${2:-30000}
ACTS=${3:-"0.050 0.055 0.060 0.065 0.070 0.075 0.080 0.085 0.090"}
cd "$(dirname "$0")"
i=0
for act in $ACTS; do
    i=$((i+1))
    d=L${L}/act_${act}
    mkdir -p $d
    cat > $d/input.json <<EOJ
{
    "domain":           [${L}, ${L}],
    "simSteps":         ${STEPS},
    "warmUp":           0,
    "seed":             $((L*100+i)),
    "collOp":           21,
    "lc":               1,
    "species":          [ { "dens": 20, "act": ${act} } ],
    "domainWalls":      1,
    "defectsOut":       20,
    "checkpointOut":    5000,
    "synopsisOut":      1
}
EOJ
done
