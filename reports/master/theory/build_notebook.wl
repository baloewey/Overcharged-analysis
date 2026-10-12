(* Builds ShieldingTheory.nb (text + live input cells), then evaluates it with the front end and saves it with outputs.
   Run: wolframscript -file build_notebook.wl *)
dir = DirectoryName[$InputFileName];
root = FileNameJoin[{dir, "..", "..", ".."}];
T[s_] := Cell[s, "Text"]; S[s_] := Cell[s, "Section"]; SS[s_] := Cell[s, "Subsection"]; I0[s_] := Cell[s, "Input"];
M[s_] := Cell[s, "Text", FontSlant -> "Italic", FontColor -> GrayLevel[0.35]];

cells = {
 Cell["Overcharged cavities: elastic theory of the wall layer (shielding)", "Title"],
 Cell["Benjamin Loewe (with Claude), 2026-10-06. Companion to Sec. \"Elastic theory of the wall layer\" of reports/master/master.pdf. Built by build_notebook.wl from the scripts shielding_disk.wl and shielding_ellipse.wl.", "Subtitle"],
 T["How to use: evaluate top to bottom (Evaluation > Evaluate Notebook). Energies are in units of \[Pi]K, lengths in MPCD cells. The ellipse relaxation (Sec. 5) takes ~10 min, so by default it loads the saved relaxed positions; set recompute = True to redo it."],
 I0["SetDirectory[" <> ToString[dir, InputForm] <> "];\nprojectRoot = " <> ToString[root, InputForm] <> ";"],

 S["1. Model"],
 T["One-constant Frank elasticity, director angle \[Theta], defects of charge q_i at z_i (complex positions). Use the field \[Psi] conjugate to \[Theta] (\[Del]\[Psi] = z\[Cross]\[Del]\[Theta]): defects are point sources, \[Del]\[Psi]^2 = 2\[Pi] \[Sum] q_i \[Delta](z - z_i). Strong anchoring fixes \[Theta] along the wall, so \[PartialD]_s\[Theta] is prescribed: a Neumann condition for \[Psi], with flux density equal to the imposed winding per length. In a disk of radius R the Neumann Green's function has images of the SAME sign at R^2/z\[Conjugate], so the wall repels every defect. Energy (core terms dropped):"],
 T["E/(\[Pi]K) = - \[Sum]_{i\[NotEqual]j} q_i q_j Log|z_i - z_j| - \[Sum]_{i,j} q_i q_j Log|R^2 - z_i Conjugate[z_j]|"],
 T["Pieces: -2\[Pi]K q_i q_j Log r_ij is -(\[Pi]/2)K Log r for like \[PlusMinus]1/2 (repulsion); the self-image term -\[Pi]K q_i^2 Log(R^2 - r_i^2) = -(\[Pi]/4)K Log(R^2 - r^2) is the repulsion from the wall. These are the poster's U_ij and U_W,i with the poster's sign convention U = -E; the formula above adds the cross-image terms."],
 I0["Edisk[zs_List, qs_List, R_] := Module[{n = Length[zs]},\n  -Sum[qs[[i]] qs[[j]] Log[Abs[zs[[i]] - zs[[j]]]], {i, n}, {j, Complement[Range[n], {i}]}]\n  - Sum[qs[[i]] qs[[j]] Log[Abs[R^2 - zs[[i]] Conjugate[zs[[j]]]]], {i, n}, {j, n}]];\n(* single defect: energy against distance from the centre -> repulsion from the wall *)\nPlot[Edisk[{r}, {1/2}, 1] - Edisk[{0.}, {1/2}, 1], {r, 0, 0.98}, Frame -> True,\n FrameLabel -> {\"r/R\", \"E (\[Pi]K)\"}, PlotLabel -> \"one +1/2 in a disk: the wall repels\"]"],

 S["2. Exact shielding in a circle"],
 T["For N wall charges q_w evenly spaced on a ring of radius r_w: \[Sum]_j Log|z - r_w e^(2\[Pi]i j/N)| = Log|z^N - r_w^N|. A bulk charge q_t at z = r e^(i\[CurlyPhi]) then feels"],
 T["V/(\[Pi]K) = const - 2 q_t q_w [ Log|1 - (r/r_w)^N e^(iN\[CurlyPhi])| + Log|1 - (r r_w/R^2)^N e^(iN\[CurlyPhi])| ] - q_t^2 Log(R^2 - r^2)"],
 I0["(* the identity, numerically *)\nWith[{z = 0.31 + 0.47 I, r = 0.9, n = 7}, Sum[Log[Abs[z - r Exp[2 Pi I j/n]]], {j, 0, n - 1}] - Log[Abs[z^n - r^n]]]"],
 T["Both brackets vanish on angular average (2D shell theorem, images included): inside the ring the bulk does not feel the wall layer, whatever its charge. Only a corrugation of relative size (r/r_w)^N remains near the wall."],
 I0["Table[NIntegrate[Log[Abs[1 - a Exp[I t]]], {t, 0, 2 Pi}]/(2 Pi), {a, {0.3, 0.7, 0.95}}] // Chop"],
 I0["Vring[r_, phi_, nn_, rw_, R_, qt_, qw_] := -2 qt qw (nn Log[rw] + Log[Abs[1 - (r/rw)^nn Exp[I nn phi]]] +\n    2 nn Log[R] + Log[Abs[1 - (r rw/R^2)^nn Exp[I nn phi]]]);\nVself[r_, R_, qt_] := -qt^2 Log[R^2 - r^2];\nManipulate[\n DensityPlot[If[x^2 + y^2 < 1, Vring[Sqrt[x^2 + y^2], Arg[x + I y], nn, rw, 1, 1/2, qw], None], {x, -1, 1}, {y, -1, 1},\n  PlotPoints -> 90, ColorFunction -> \"TemperatureMap\", PlotRange -> {-1, 1}, ClippingStyle -> Automatic, PlotLegends -> Automatic,\n  PlotLabel -> \"energy of a bulk +1/2 from the wall ring (\[Pi]K)\"],\n {{nn, 14, \"wall defects N\"}, 4, 40, 2}, {{rw, 0.94, \"ring radius r_w/R\"}, 0.8, 0.99}, {{qw, -1/2, \"wall charge\"}, {-1/2, 1/2}}]"],

 S["3. Where the passive wall layer sits"],
 T["Energy of a passive ring of N charges q at radius r: E/(\[Pi]K) = -q^2 N [(N-1) Log r + Log(R^(2N) - r^(2N))] + const. Claimed minimum r* = R ((N-1)/(3N-1))^(1/(2N)); for large N the depth is (Log 3/(4\[Pi])) \[Times] spacing \[TildeTilde] 0.087 \[Times] spacing."],
 I0["Ering[r_, nn_, R_, q_] := -q^2 nn ((nn - 1) Log[r] + Log[R^(2 nn) - r^(2 nn)]);\nrstar = R ((n - 1)/(3 n - 1))^(1/(2 n));\nFullSimplify[D[Ering[r, n, R, q], r] /. r -> rstar, n > 1 && R > 0 && q > 0]   (* = 0: stationary *)"],
 I0["D[Ering[r, 12, 1, 1/2], {r, 2}] /. r -> N[(11/35)^(1/24)]   (* > 0: minimum *)"],
 I0["(* with a bulk charge Qb at the centre (active layer) *)\nEact[r_, nn_, R_, q_, Qb_] := Ering[r, nn, R, q] - 2 q Qb nn Log[r];\nFullSimplify[D[Eact[r, n, R, q, Qb], r] /. r -> R ((n - 1 + 2 Qb/q)/(3 n - 1 + 2 Qb/q))^(1/(2 n)), n > 1 && R > 0 && q > 0 && Qb >= 0]"],
 I0["Plot[(1 - ((nn - 1)/(3 nn - 1))^(1/(2 nn))) nn/(2 Pi), {nn, 2, 60}, Frame -> True, FrameLabel -> {\"N\", \"depth / spacing\"},\n Epilog -> {Dashed, Line[{{0, Log[3]/(4 Pi)}, {60, Log[3]/(4 Pi)}}]}, PlotLabel -> \"depth in units of the spacing (dashed: Log 3/4\[Pi])\"]"],
 SS["Comparison with all passive circles (no free parameter)"],
 T["Data: reports/master/data/passive_circles.tsv (468 passive circles of the main sweep; s_mean = mean distance of the imposed defects from the wall). The fit gives slope \[TildeTilde] 0.98 and a constant offset of \[TildeTilde] 0.6 cells (where the cell-based anchoring puts the effective wall)."],
 I0["data = Import[FileNameJoin[{projectRoot, \"reports\", \"master\", \"data\", \"passive_circles.tsv\"}], \"TSV\"];\nrows = Select[Rest[data], #[[2]] <= 0.151 && Abs[#[[3]]] >= 1 &];\npred = (#[[4]] (1 - ((2 Abs[#[[3]]] - 1)/(6 Abs[#[[3]]] - 1))^(1/(4 Abs[#[[3]]])))) & /@ rows;\nmeas = rows[[All, 8]];\nfit = LinearModelFit[Transpose[{pred, meas}], x, x];\n{Length[rows], fit[\"BestFitParameters\"], Correlation[pred, meas]}"],
 I0["Show[ListPlot[Transpose[{pred, meas}], Frame -> True, FrameLabel -> {\"predicted depth (cells)\", \"measured depth (cells)\"}],\n Plot[{x, fit[x]}, {x, 0, 3.6}, PlotStyle -> {Dashed, Black}]]"],

 S["4. A bulk defect approaching the wall layer"],
 T["Mirror circle (k = -6, R = 31.83 = 1.43 \[ScriptL]_a): 2|k|+2 = 14 wall -1/2 at the predicted depth plus the 0.61-cell offset. A bulk +1/2 moves along a radius towards a wall defect, or towards the gap between two; compared with a free -1/2."],
 I0["approach[R0_, nW_] := Module[{rw = R0 - (0.61 + R0 (1 - ((nW - 1)/(3 nW - 1))^(1/(2 nW))))},\n  {Function[d, Vring[rw - d, 0, nW, rw, R0, 1/2, -1/2] + Vself[rw - d, R0, 1/2]],\n   Function[d, Vring[rw - d, Pi/nW, nW, rw, R0, 1/2, -1/2] + Vself[rw - d, R0, 1/2]], rw}];\nVfree[d_] := -2 (1/2) (-1/2) Log[d];\n{Vline, Vmid, rw} = approach[31.83, 14];\nGrid[{{\"\", \"free -1/2\", \"towards a wall -1/2\"},\n  {\"energy released 8 -> 4 cells\", Vfree[8] - Vfree[4], Vline[8] - Vline[4]},\n  {\"energy released 4 -> 1 cell\", Vfree[4] - Vfree[1], Vline[4] - Vline[1]}}, Frame -> All]"],
 I0["Manipulate[Module[{vl, vm, r},\n  {vl, vm, r} = approach[R0, nW];\n  Plot[{Vfree[d] - Vfree[8], vl[d] - vl[8], vm[d] - vm[8]}, {d, 0.8, 12}, Frame -> True, PlotRange -> {-0.8, 0.8},\n   FrameLabel -> {\"distance from the wall layer (cells)\", \"energy (\[Pi]K), zero at 8\"},\n   PlotLegends -> {\"free -1/2\", \"towards a wall -1/2\", \"between two wall -1/2\"}]],\n {{R0, 31.83, \"R (cells)\"}, 15, 64}, {{nW, 14, \"wall sites 2|k|+2\"}, 6, 40, 2}]"],

 S["5. Ellipses: conformal map"],
 T["Ellipse x^2/A^2 + y^2/B^2 < 1, c = Sqrt[A^2 - B^2]. The map w(z) = m^(1/4) JacobiSN[(2K(m)/\[Pi]) ArcSin[z/c], m], nome ((A-B)/(A+B))^2, sends it to the unit disk. With winding uniform in arclength (eMode 2) the Neumann Green's function is G(z,\[Zeta]) = Log|w-\[Omega]| + Log|1 - w \[Omega]\[Conjugate]| + H(z) + H(\[Zeta]), where H is harmonic with \[PartialD]_n H = 2\[Pi]/P - |w'(z)| on the wall (a smooth datum, solved by a Fourier series on the disk). Energy:"],
 T["E/(\[Pi]K) = -\[Sum]_{i<j} 2 q_i q_j (Log|w_i - w_j| + Log|1 - w_i w_j\[Conjugate]|) - \[Sum]_i q_i^2 (Log|w'(z_i)| + Log(1 - |w_i|^2)) - 2 Q \[Sum]_i q_i H(z_i)"],
 I0["setup[A_, B_] := Module[{c, q, m, K, wf, wp, zf, P, nF = 256, al, zs, g, cf},\n  c = Sqrt[A^2 - B^2]; q = ((A - B)/(A + B))^2; m = InverseEllipticNomeQ[q]; K = EllipticK[m];\n  wf = Function[z, m^(1/4) JacobiSN[(2 K/Pi) ArcSin[z/c], m]];\n  wp = Function[z, Evaluate[D[m^(1/4) JacobiSN[(2 K/Pi) ArcSin[zz/c], m], zz] /. zz -> z]];\n  zf = Function[w, c Sin[(Pi/(2 K)) InverseJacobiSN[w/m^(1/4), m]]];\n  P = 4 A EllipticE[1 - B^2/A^2];\n  al = Table[2 Pi j/nF, {j, 0, nF - 1}]; zs = zf /@ Exp[I al];\n  g = (2 Pi/P)/Abs[wp /@ zs] - 1; cf = Fourier[g, FourierParameters -> {-1, -1}];\n  <|\"A\" -> A, \"B\" -> B, \"w\" -> wf, \"wp\" -> wp, \"z\" -> zf, \"P\" -> P,\n    \"H\" -> Function[w, Re[Sum[2 cf[[n + 1]] w^n/n, {n, 1, nF/2 - 1}]]]|>];\nenergyE[d_, zs_, qs_] := Module[{ws = d[\"w\"] /@ zs, Q = Total[qs], n = Length[zs]},\n  -Sum[2 qs[[i]] qs[[j]] (Log[Abs[ws[[i]] - ws[[j]]]] + Log[Abs[1 - ws[[i]] Conjugate[ws[[j]]]]]), {i, 1, n - 1}, {j, i + 1, n}]\n  - Sum[qs[[i]]^2 (Log[Abs[d[\"wp\"][zs[[i]]]]] + Log[1 - Abs[ws[[i]]]^2]), {i, n}] - 2 Q Sum[qs[[i]] d[\"H\"][ws[[i]]], {i, n}]];\n(* check: the wall maps onto the unit circle *)\nWith[{d = setup[40., 20.]}, Max[Abs[Abs[d[\"w\"] /@ Table[40 Cos[t] + 20 I Sin[t], {t, 0.1, 6.2, 0.5}]] - 1]]]"],
 SS["Passive equilibrium of 12 +1/2 (k = 6), area 2500"],
 I0["recompute = False;   (* True: relax again (~10 min); False: load relaxed positions saved by shielding_ellipse.wl *)\nrelax[d_, nn_, q_, ph_] := Module[{vars, start, obj, sol},\n  vars = Flatten[Table[{u[i], p[i]}, {i, nn}]];\n  start = Flatten[Table[{{u[i], 3.0}, {p[i], 2 Pi (i - ph)/nn}}, {i, nn}], 1];\n  obj[v__] := Module[{vv = {v}}, energyE[d, d[\"z\"] /@ Table[Tanh[vv[[2 i - 1]]]^(1/2) Exp[I vv[[2 i]]], {i, nn}], ConstantArray[q, nn]]];\n  sol = FindMinimum[obj @@ vars, start, Method -> \"QuasiNewton\", MaxIterations -> 4000];\n  {Table[d[\"z\"][Tanh[u[i]]^(1/2) Exp[I p[i]]] /. sol[[2]], {i, nn}], sol[[1]]}];\nres = Association[Table[eps -> Module[{A = Sqrt[2500. eps/Pi], d},\n     d = setup[A, A/eps];\n     <|\"d\" -> d, \"zs\" -> If[recompute, First[MinimalBy[{relax[d, 12, 1/2, 0.5], relax[d, 12, 1/2, 0.]}, Last]][[1]],\n         (#[[1]] + I #[[2]]) & /@ Import[\"relaxed_eps\" <> ToString[NumberForm[If[eps == 1.0001, 1., eps], {3, 2}]] <> \".csv\"]]|>],\n   {eps, {1.0001, 1.2, 1.5, 2.0}}]];\nGraphicsRow[Table[With[{d = res[e][\"d\"], zs = res[e][\"zs\"]},\n   Graphics[{Circle[{0, 0}, {d[\"A\"], d[\"B\"]}], Red, PointSize[0.03], Point[ReIm /@ zs]}, PlotLabel -> \"\[Epsilon] = \" <> ToString[e]]],\n  {e, Keys[res]}], ImageSize -> 800]"],
 T["The passive layer stays evenly spaced in arclength (gaps at the tips = gaps on the flanks), as in the passive simulations; it does not adopt the screening (conductor) distribution, which would need gaps \[Epsilon] times smaller at the tips."],
 SS["Leftover field inside: is the ellipse screened?"],
 I0["layerE[d_, zs_, qs_, z_, qt_] := Module[{w = d[\"w\"][z], ws = d[\"w\"] /@ zs},\n  -2 qt Sum[qs[[i]] (Log[Abs[w - ws[[i]]]] + Log[Abs[1 - w Conjugate[ws[[i]]]]]), {i, Length[zs]}] - 2 Total[qs] qt d[\"H\"][w]];\nTable[Module[{d = res[e][\"d\"], zs = res[e][\"zs\"], A, B, pts, lp},\n   A = d[\"A\"]; B = d[\"B\"];\n   pts = Select[Flatten[Table[x + I y, {x, -A, A, A/14}, {y, -B, B, B/10}]], (Re[#]/A)^2 + (Im[#]/B)^2 < 0.6^2 &];\n   lp = layerE[d, zs, ConstantArray[1/2, 12], #, 1/2] & /@ pts;\n   {e, Max[lp] - Min[lp]}], {e, Keys[res]}] // TableForm[#, TableHeadings -> {None, {\"\[Epsilon]\", \"layer energy range inside the 0.6 level (\[Pi]K)\"}}] &"],
 T["The range grows ~20\[Times] from the circle to \[Epsilon] = 2 but stays \[LessTilde] 0.06 \[Pi]K in the inner 60% (Sec. 6: \[Pi]K \[TildeTilde] 100 k_BT, so \[TildeTilde] 6 k_BT): the shielding is imperfect in ellipses but modest compared with the approach energies near the wall."],

 S["6. Energy scale"],
 T["Measured from equilibrium director fluctuations of passive defect-free periodic boxes (runs/Kcalib, frank_constant.py; density 20, mean-field potential 10 k_BT): K1 (splay) \[TildeTilde] 42 k_BT, K3 (bend) \[TildeTilde] 23 k_BT, independent of box size (64^2 and 96^2). One-constant K \[TildeTilde] 32 k_BT, \[Pi]K \[TildeTilde] 100 k_BT."],
 I0["piK = Pi 32.4;   (* k_BT *)\n{\"barrier towards a wall -1/2\" -> 0.04 piK, \"pull 4 -> 1 cell (wall / free)\" -> {0.30, 0.69} piK,\n \"leftover field inside, eps = 2, inner 60%\" -> 0.059 piK, \"Sqrt[K/zeta] (cells), zeta = 0.07\" -> Sqrt[32.4/0.07]}"],

 S["7. Limits"],
 T["One-constant elasticity (the measured 2D constants differ by a factor ~1.8); orientation-dependent interactions ignored; core energies and the cell scale enter only through the empirical 0.61-cell offset; no activity or flow. The negative-winding ellipses (master document) show that the strong loss of the layer with positive winding needs self-propelled wall defects; the elastic shielding weakens only modestly with elongation."]
};

nbExpr = Notebook[cells, WindowTitle -> "Shielding theory"];
file = FileNameJoin[{dir, "ShieldingTheory.nb"}];
Export[file, nbExpr];
Print["unevaluated notebook written: ", file];
ok = Quiet@Check[UsingFrontEnd[
   Module[{nb = NotebookOpen[file, Visible -> False]},
    NotebookEvaluate[nb, InsertResults -> True];
    NotebookSave[nb]; NotebookClose[nb]]; True], False];
Print["evaluated with the front end: ", ok];
