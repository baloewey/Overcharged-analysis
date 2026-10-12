#!/bin/bash
# Passive (act=0) validation of the customAnch hack in a circular cavity.
# R=15, box 33x33 (smallest odd L >= 2R+2), centre on a cell centre.
# Runs: k1_ca0 (customAnch off -> true normal), k1, k2, k3, k5.
cd "$(dirname "$0")"
R=15; L=33; C=16.5; POP=14137      # POP = round(20*pi*R^2)
i=0
for name in k1_ca0 k1 k2 k3 k5; do
    i=$((i+1))
    ca=${name#k}; ca=${ca%_ca0}; [ "$name" = k1_ca0 ] && ca=0
    mkdir -p $name
    cat > $name/input.json <<EOJ
{
    "domain":           [${L}, ${L}],
    "simSteps":         20000,
    "warmUp":           0,
    "seed":             $((500+i)),
    "collOp":           21,
    "lc":               1,
    "species":          [ { "pop": ${POP}, "act": 0.0 } ],
    "BC":
    [
        {
            "comment":      "No-slip circular cavity, homeotropic custom anchoring",
            "phantom":      1,
            "Q":            [${C}, ${C}, 0],
            "aInv":         [1, 1, 1],
            "P":            [2, 2, 2, 2],
            "R":            ${R},
            "DN":           0,
            "MVN":          -1,
            "MVT":          -1,
            "MUN":          1,
            "MUT":          0,
            "inv":          1,
            "customAnch":   ${ca}
        }
    ],
    "defectsOut":       100,
    "dirSOut":          1000,
    "topoFieldOut":     1000,
    "synopsisOut":      1
}
EOJ
done
