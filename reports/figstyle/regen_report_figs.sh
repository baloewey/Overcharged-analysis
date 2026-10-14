#!/bin/bash
# Regenerate every figure of the theory report with the house style (reports/figstyle/sitecustomize.py; PNG + vector PDF) and copy them into theory_report/figs
# and presentation/figs (where used). Usage: reports/figstyle/regen_report_figs.sh   (log: reports/figstyle/regen.log)
O=/home/baloewe/Physics/Overcharged; T=$O/reports/master/theory; E=$T/energetics; F=$O/reports/theory_report/figs; P=$O/reports/presentation/figs
export PYTHONPATH=$O/reports/figstyle; export NPROC=4
r() { local d=$1; shift; echo "$(date +%T) $*"; (cd "$d" && nice -n 10 "$@" > /dev/null 2>> $O/reports/figstyle/regen.err) || echo "   FAILED: $*"; }
r $E python3 plot_etsg.py
r $E python3 plot_etsg_reg.py
r $E env REPLOT=1 python3 perlength.py
r $E env REPLOT=1 python3 vortdef.py ff36
r $E env REPLOT=1 SIGMA=2.5 python3 qcrit.py ff36s25
r $E env REPLOT=1 python3 jetmap.py ff12
r $E env REPLOT=1 TAG=_ff12 python3 selfflow.py   # from selfflow_ff12.txt; full recompute: TAG=_ff12 selfflow.py <the 12 fineflow/fineflow2 runs>
r $E python3 compare_bulk.py
r $O/runs python3 $E/tiltlaw.py $(ls -d $O/runs/lambda_circle/*/)
r $T/dirstokes python3 fig_dir.py
r $T/effbc python3 effective_bc.py
r $T/ejection python3 tilted_structures.py
r $T/shear_tilt python3 fig.py
r $T/shapes python3 phase_diagram.py
r $T/selfprop python3 fig_restoring.py
r $T/channel python3 channel_fig.py
r $T/channel python3 layer_profile.py
r $O/reports/theory_report python3 fig_flatwall.py
r $O/reports/theory_report python3 fig_cont_ellipse.py
r $O/reports/theory_report python3 make_figs.py
r $O/reports/master/scripts python3 fig_distributions.py --workers 4
for f in $E/fig_etsg_shapes.png $E/fig_etsg_reg.png $E/fig_perlength.png $E/fig_vortdef_ff36.png $E/fig_qcrit_ff36s25.png $E/fig_jetmap_ff12.png $E/fig_selfflow_ff12.png \
         $E/fig_compare_bulk.png $E/fig_tiltlaw.png $T/dirstokes/dirstokes.png $T/effbc/effective_bc.png $T/ejection/tilted_structures.png $T/shear_tilt/shear_tilt.png \
         $T/shapes/phase_shapes.png $T/channel/fig_channel.png $T/channel/fig_layer_profile.png $O/reports/master/figs/maps_circles.png; do
    cp $f $F/; [ -f ${f%.png}.pdf ] && cp ${f%.png}.pdf $F/; b=$(basename $f); [ -f $P/$b ] && cp $f $P/$b
done
for b in fig_restoring.png fig_flatwall.png fig_cont_ellipse.png fig_lock.png fig_ejection.png fig_shapes_mech.png; do [ -f $P/$b ] && cp $F/$b $P/$b; done
echo "$(date +%T) done"
