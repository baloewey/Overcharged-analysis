#!/bin/bash
# One-constant K = 44 reruns of the escape/steady-state energetics (2026-10-22; binding.Ellipse reads FRANK_K1/K3). Office, nice 19, after the user's jobs as they free cores.
cd /home/baloewe/Physics/Overcharged/reports/master/theory/energetics
export FRANK_K1=44 FRANK_K3=44 NPROC=6
R=/home/baloewe/Physics/Overcharged/runs; SW=$R/shapes_server/oval
nice -n 19 python3 steadyE.py A_K44 $R/fineflow/eps1.0_s1 $R/fineflow2/eps1.0_s3 $R/fineflow2/eps1.4_s3 $R/fineflow2/eps1.4_s4 $R/fineflow/eps1.6_s1 $R/fineflow2/eps1.6_s3 $R/fineflow/eps2.0_s1 $R/fineflow2/eps2.0_s3 $R/lambda/lam4.0_eps2.0_s1 $R/lambda/lam4.0_eps2.0_s2 $R/lambda/lam1.2_eps2.0_s1 $R/lambda/lam1.2_eps2.0_s2 > steadyE_A_K44.log 2>&1
nice -n 19 python3 steadyE.py oval_K44 $R/ovalfine/ba* $SW/ba2.50/Req1.5/k06/s0* $SW/ba1.80/Req1.5/k06/s0* > steadyE_oval_K44.log 2>&1
SW2=$R/shapes_server
nice -n 19 python3 energetics_general.py shapes_reg_K44 $(ls -d $SW2/ellipse/*/*/*/s0* $SW2/circle/circle/*/*/s0*) > etsg_shapes_reg_K44.log 2>&1
echo "$(date) done" >> queue_K44.log
