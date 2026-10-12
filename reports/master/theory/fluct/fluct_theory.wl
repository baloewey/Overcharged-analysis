(* Stage 1: harmonic (Hessian) theory of the passive wall-lattice fluctuations, no free parameter.
   Energy in units of pi K (Sec. theory): disk  E = - Sum_{i!=j} qi qj Log|zi - zj| - Sum_{i,j} qi qj Log|R^2 - zi Conj[zj]|;
   ellipse: conformal map + harmonic correction (shielding_ellipse.wl). q = 1/2, N = 2k passive +1/2.
   Equipartition: <dx dx^T> = (kBT/(pi K)) H^+, H = Hessian of E in Cartesian coordinates, pseudo-inverse without the rotation
   zero mode of the disk. K = 32.5 kBT (mean of K1 = 42.0, K3 = 23.1; Sec. frank) -> kBT/(pi K) = 0.0098.
   Output: rms tangential / normal displacement per defect, spectrum per Fourier mode (disk), lowest mode (ellipse: the sliding mode).
   Run: wolframscript -file fluct_theory.wl -> fluct_theory.json *)
dir = DirectoryName[$InputFileName];
kT = 1/(Pi 32.5); q = 1/2;
(* ---------- disk ---------- *)
diskTheory[R_, nn_] := Module[{rs = R ((nn - 1)/(3 nn - 1))^(1/(2 nn)), ph, xs, vars, E, H, es, vs, C, tang, norm, rule, um, specT, specN},
  ph = Table[2 Pi (j - 1/2)/nn, {j, nn}];
  vars = Flatten[Table[{X[j], Y[j]}, {j, nn}]];
  E = (-Sum[If[i != j, q^2 Log[(X[i] - X[j])^2 + (Y[i] - Y[j])^2]/2, 0], {i, nn}, {j, nn}]
       - Sum[q^2 Log[(R^2 - X[i] X[j] - Y[i] Y[j])^2 + (X[i] Y[j] - Y[i] X[j])^2]/2, {i, nn}, {j, nn}]);
  rule = Flatten[Table[{X[j] -> rs Cos[ph[[j]]], Y[j] -> rs Sin[ph[[j]]]}, {j, nn}]];
  H = N[D[E, {vars, 2}] /. rule];
  {es, vs} = Eigensystem[H];
  C = kT Total[Table[If[Abs[es[[m]]] > 10^-8 Max[Abs[es]], Outer[Times, vs[[m]], vs[[m]]]/es[[m]], 0 H], {m, Length[es]}]];
  tang = Table[With[{t = ph[[j]]}, ReplacePart[ConstantArray[0., 2 nn], {2 j - 1 -> -Sin[t], 2 j -> Cos[t]}]], {j, nn}];
  norm = Table[With[{t = ph[[j]]}, ReplacePart[ConstantArray[0., 2 nn], {2 j - 1 -> -Cos[t], 2 j -> -Sin[t]}]], {j, nn}];   (* towards the wall... sign irrelevant *)
  specT = Table[With[{v = Sum[Exp[-2 Pi I m (j - 1)/nn] tang[[j]], {j, nn}]/Sqrt[nn]}, Re[Conjugate[v] . C . v]], {m, 0, Floor[nn/2]}];
  specN = Table[With[{v = Sum[Exp[-2 Pi I m (j - 1)/nn] norm[[j]], {j, nn}]/Sqrt[nn]}, Re[Conjugate[v] . C . v]], {m, 0, Floor[nn/2]}];
  <|"N" -> nn, "R" -> R, "depth" -> R - rs, "sig_t" -> Sqrt[Mean[Table[tang[[j]] . C . tang[[j]], {j, nn}]]],
    "sig_n" -> Sqrt[Mean[Table[norm[[j]] . C . norm[[j]], {j, nn}]]], "spec_t" -> specT, "spec_n" -> specN,
    "lowest_nonzero" -> Min[Select[es, # > 10^-8 Max[Abs[es]] &]]|>];
res = <||>;
Do[res["circle_R31.83_k" <> IntegerString[k, 10, 2]] = diskTheory[31.830989, 2 k]; Print["disk k ", k, ": ", KeyDrop[res["circle_R31.83_k" <> IntegerString[k, 10, 2]], {"spec_t", "spec_n"}]], {k, 2, 10}];
res["circle_R27.28_k06"] = diskTheory[27.283705, 12]; Print["disk R27.28: ", KeyDrop[res["circle_R27.28_k06"], {"spec_t", "spec_n"}]];
(* ---------- ellipse (functions from shielding_ellipse.wl) ---------- *)
setup[A_, B_] := Module[{c, qq, m, K, wf, wp, zf, P, nF = 256, al, zs, g, cf, Hfun},
  c = Sqrt[A^2 - B^2]; qq = ((A - B)/(A + B))^2; m = InverseEllipticNomeQ[qq]; K = EllipticK[m];
  wf = Function[z, m^(1/4) JacobiSN[(2 K/Pi) ArcSin[z/c], m]];
  wp = Function[z, Evaluate[D[m^(1/4) JacobiSN[(2 K/Pi) ArcSin[zz/c], m], zz] /. zz -> z]];
  zf = Function[w, c Sin[(Pi/(2 K)) InverseJacobiSN[w/m^(1/4), m]]];
  P = 4 A EllipticE[1 - B^2/A^2];
  al = Table[2 Pi j/nF, {j, 0, nF - 1}]; zs = zf /@ Exp[I al];
  g = (2 Pi/P)/Abs[wp /@ zs] - 1; cf = Fourier[g, FourierParameters -> {-1, -1}];
  Hfun = Function[w, Re[Sum[2 cf[[n + 1]] w^n/n, {n, 1, nF/2 - 1}]]];
  <|"A" -> A, "B" -> B, "w" -> wf, "wp" -> wp, "z" -> zf, "P" -> P, "H" -> Hfun|>];
energy[d_, zs_, qs_] := Module[{ws = d["w"] /@ zs, Q = Total[qs], n = Length[zs]},
  -Sum[2 qs[[i]] qs[[j]] (Log[Abs[ws[[i]] - ws[[j]]]] + Log[Abs[1 - ws[[i]] Conjugate[ws[[j]]]]]), {i, 1, n - 1}, {j, i + 1, n}]
  - Sum[qs[[i]]^2 (Log[Abs[d["wp"][zs[[i]]]]] + Log[1 - Abs[ws[[i]]]^2]), {i, n}] - 2 Q Sum[qs[[i]] d["H"][ws[[i]]], {i, n}]];
ellTheory[eps_, tag_] := Module[{A = Sqrt[2500. eps/Pi], B, d, z0, nn, f, x0, sol, x1, h = 0.004, H, es, vs, C, Cint, tang, norm, ph, tt, low, tw},
  B = A/eps; d = setup[A, B];
  z0 = (#[[1]] + I #[[2]]) & /@ Import[FileNameJoin[{ParentDirectory[dir], "relaxed_eps" <> tag <> ".csv"}]]; nn = Length[z0];
  f[x_?(VectorQ[#, NumericQ] &)] := energy[d, Table[x[[2 j - 1]] + I x[[2 j]], {j, nn}], ConstantArray[q, nn]];
  x0 = Flatten[ReIm /@ z0];
  sol = FindMinimum[f[Array[xx, 2 nn]], Transpose[{Array[xx, 2 nn], x0}], Method -> "QuasiNewton", AccuracyGoal -> 10, PrecisionGoal -> 10, MaxIterations -> 2000];
  x1 = Array[xx, 2 nn] /. sol[[2]];
  H = Table[If[a <= b, (f[x1 + h (UnitVector[2 nn, a] + UnitVector[2 nn, b])] - f[x1 + h (UnitVector[2 nn, a] - UnitVector[2 nn, b])]
        - f[x1 - h (UnitVector[2 nn, a] - UnitVector[2 nn, b])] + f[x1 - h (UnitVector[2 nn, a] + UnitVector[2 nn, b])])/(4 h^2), 0], {a, 2 nn}, {b, 2 nn}];
  H = H + Transpose[H] - DiagonalMatrix[Diagonal[H]];
  {es, vs} = Eigensystem[H]; C = kT Inverse[H];
  (* local wall tangent/normal at the nearest wall point *)
  tt = Table[With[{z = x1[[2 j - 1]] + I x1[[2 j]]}, tw /. Last[NMinimize[{Abs[z - (A Cos[tw] + I B Sin[tw])]^2, -Pi <= tw <= Pi}, tw]]], {j, nn}];
  tang = Table[With[{t = tt[[j]], v = Normalize[{-A Sin[tt[[j]]], B Cos[tt[[j]]]}]}, ReplacePart[ConstantArray[0., 2 nn], {2 j - 1 -> v[[1]], 2 j -> v[[2]]}]], {j, nn}];
  norm = Table[With[{v = Normalize[{B Cos[tt[[j]]], A Sin[tt[[j]]]}]}, ReplacePart[ConstantArray[0., 2 nn], {2 j - 1 -> v[[1]], 2 j -> v[[2]]}]], {j, nn}];
  low = Last[Ordering[es, 1]];
  Cint = kT Total[Table[If[m == low, 0 H, Outer[Times, vs[[m]], vs[[m]]]/es[[m]]], {m, 2 nn}]];   (* without the sliding mode *)
  <|"N" -> nn, "eps" -> eps, "x" -> Partition[x1, 2], "relax_shift" -> Max[Abs[x1 - x0]], "min_eig" -> Min[es],
    "sig_t_each" -> Sqrt[Table[tang[[j]] . C . tang[[j]], {j, nn}]], "sig_n_each" -> Sqrt[Table[norm[[j]] . C . norm[[j]], {j, nn}]],
    "sig_t" -> Sqrt[Mean[Table[tang[[j]] . C . tang[[j]], {j, nn}]]], "sig_n" -> Sqrt[Mean[Table[norm[[j]] . C . norm[[j]], {j, nn}]]],
    "sig_coll" -> Sqrt[(Total[tang]/nn) . C . (Total[tang]/nn)],
    "sig_t_int_each" -> Sqrt[Table[tang[[j]] . Cint . tang[[j]], {j, nn}]], "sig_n_int_each" -> Sqrt[Table[norm[[j]] . Cint . norm[[j]], {j, nn}]],
    "sig_t_int" -> Sqrt[Mean[Table[tang[[j]] . Cint . tang[[j]], {j, nn}]]], "sig_n_int" -> Sqrt[Mean[Table[norm[[j]] . Cint . norm[[j]], {j, nn}]]],
    "wall_t" -> tt, "second_eig" -> Sort[es][[2]],
    "sliding_eig" -> Min[es], "sliding_overlap" -> Abs[Normalize[Total[tang]] . vs[[Last[Ordering[es, 1]]]]], "cov" -> C|>];
Do[With[{r = ellTheory[ToExpression[tg], tg]}, res["ellipse_eps" <> StringDrop[tg, -1] <> "_k06"] = r;
   Print["ellipse ", tg, ": ", KeyDrop[r, {"x", "cov", "sig_t_each", "sig_n_each", "sig_t_int_each", "sig_n_int_each", "wall_t"}], " internal tang each ", NumberForm[r["sig_t_int_each"], 3]]],
  {tg, {"1.20", "1.50", "2.00"}}];
Export[FileNameJoin[{dir, "fluct_theory.json"}], res, "RawJSON"];
