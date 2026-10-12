(* Defect electrostatics in an ellipse with strong anchoring and winding uniform in arclength (eMode 2), via conformal mapping.
   Ellipse x^2/A^2 + y^2/B^2 < 1, c = Sqrt[A^2 - B^2]. Map to the unit disk: w(z) = m^(1/4) JacobiSN[(2 K(m)/Pi) ArcSin[z/c], m],
   nome q = ((A - B)/(A + B))^2, m = InverseEllipticNomeQ[q].
   Neumann Green's function with uniform flux 2 Pi/P on the ellipse:
     G(z, zeta) = Log|w - omega| + Log|1 - w Conjugate[omega]| + H(z) + H(zeta),
   H harmonic with dH/dn = 2 Pi/P - |w'(z)| on the wall (the disk Green's function carries flux |w'| per ellipse length).
   On the disk: dHt/dr = (2 Pi/P)/|w'(z(alpha))| - 1 at r = 1 -> Fourier series, Ht = Re Sum_n (2 cn/n) w^n.
   Energy (units of pi K):  E = - Sum_{i<j} 2 qi qj (Log|wi - wj| + Log|1 - wi Conj[wj]|)
                                - Sum_i qi^2 (Log|w'(zi)| + Log[1 - |wi|^2]) - 2 Q Sum_i qi H(zi),   Q = Sum qi.
   Run: wolframscript -file shielding_ellipse.wl *)

dir = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{dir, "shielding_ellipse_results.txt"}]];
say[x__] := (WriteString[out, StringJoin[ToString /@ {x}], "\n"]; Print[x]);

setup[A_, B_] := Module[{c, q, m, K, wf, wp, zf, P, nF = 256, al, zs, g, cf, Hfun},
  c = Sqrt[A^2 - B^2]; q = ((A - B)/(A + B))^2; m = InverseEllipticNomeQ[q]; K = EllipticK[m];
  wf = Function[z, m^(1/4) JacobiSN[(2 K/Pi) ArcSin[z/c], m]];
  wp = Function[z, Evaluate[D[m^(1/4) JacobiSN[(2 K/Pi) ArcSin[zz/c], m], zz] /. zz -> z]];
  zf = Function[w, c Sin[(Pi/(2 K)) InverseJacobiSN[w/m^(1/4), m]]];
  P = 4 A EllipticE[1 - B^2/A^2];
  al = Table[2 Pi j/nF, {j, 0, nF - 1}];
  zs = zf /@ Exp[I al];
  g = (2 Pi/P)/Abs[wp /@ zs] - 1;                       (* radial derivative of Ht on the unit circle *)
  cf = Fourier[g, FourierParameters -> {-1, -1}];        (* g(alpha) = Sum_n cf[n] e^{i n alpha} *)
  (* Ht(w) = Re Sum_{n>=1} 2 cf[n] w^n / n  (g real: c_{-n} = conj c_n); n up to nF/2 - 1 *)
  Hfun = Function[w, Re[Sum[2 cf[[n + 1]] w^n/n, {n, 1, nF/2 - 1}]]];
  <|"A" -> A, "B" -> B, "w" -> wf, "wp" -> wp, "z" -> zf, "P" -> P, "H" -> Hfun, "g0" -> cf[[1]], "m" -> m|>];

energy[d_, zs_, qs_] := Module[{ws = d["w"] /@ zs, Q = Total[qs], n = Length[zs]},
  -Sum[2 qs[[i]] qs[[j]] (Log[Abs[ws[[i]] - ws[[j]]]] + Log[Abs[1 - ws[[i]] Conjugate[ws[[j]]]]]), {i, 1, n - 1}, {j, i + 1, n}]
  - Sum[qs[[i]]^2 (Log[Abs[d["wp"][zs[[i]]]]] + Log[1 - Abs[ws[[i]]]^2]), {i, n}]
  - 2 Q Sum[qs[[i]] d["H"][ws[[i]]], {i, n}]];

(* passive ring relaxation: positions parametrised in the disk (w = rho e^{i phi}) so the wall is never crossed *)
relax[d_, nn_, q_, ph_:0.5] := Module[{vars, start, obj, sol, zs},
  vars = Flatten[Table[{Subscript[u, i], Subscript[p, i]}, {i, nn}]];
  start = Flatten[Table[{{Subscript[u, i], 3.0}, {Subscript[p, i], 2 Pi (i - ph)/nn}}, {i, nn}], 1];
  obj[v__] := Module[{vv = {v}, ws}, ws = Table[Tanh[vv[[2 i - 1]]]^(1/2) Exp[I vv[[2 i]]], {i, nn}];   (* |w| < 1 *)
     energy[d, d["z"] /@ ws, ConstantArray[q, nn]]];
  sol = FindMinimum[obj @@ vars, start, Method -> "QuasiNewton", MaxIterations -> 4000, AccuracyGoal -> 8, PrecisionGoal -> 8];
  zs = Table[d["z"][Tanh[Subscript[u, i]]^(1/2) Exp[I Subscript[p, i]]] /. sol[[2]], {i, nn}];
  {zs, sol[[1]]}];

distToWall[d_, z_] := Module[{t}, NMinValue[{Abs[z - (d["A"] Cos[t] + I d["B"] Sin[t])], 0 <= t <= 2 Pi}, t]];
tipRatio[d_, zs_] := Module[{A = d["A"], B = d["B"], tt, W, seg, tip, near, nt, nf},
  tt = Table[2 Pi j/4000, {j, 0, 3999}]; W = A Cos[tt] + I B Sin[tt]; seg = Abs[RotateLeft[W] - W]; tip = Abs[Re[W]] > A/Sqrt[2];
  near = (First[Ordering[Abs[W - #], 1]] &) /@ zs;
  nt = Count[tip[[near]], True]; nf = Length[zs] - nt;
  (nt/Total[Pick[seg, tip]])/(nf/Total[Pick[seg, tip, False]])];

area = 2500.; nn = 12; q = 1/2;
say["Ellipses of area ", area, ", N = ", nn, " passive +1/2 (k = 6), winding uniform in arclength"];
res = <||>;
Do[
  Module[{A = Sqrt[area eps/Pi], B, d, zs, e, depths, chk},
    B = A/eps; d = setup[A, B];
    chk = Max[Abs[Abs[d["w"] /@ Table[A Cos[t] + I B Sin[t], {t, 0.1, 6.2, 0.7}]] - 1]];
    Module[{r1 = relax[d, nn, q, 0.5], r2 = relax[d, nn, q, 0.0]},
      say["   eps ", eps, ": relaxed energies for the two start orientations ", r1[[2]], ", ", r2[[2]]];
      {zs, e} = If[r1[[2]] <= r2[[2]], r1, r2]];
    depths = distToWall[d, #] & /@ zs;
    res[eps] = <|"d" -> d, "zs" -> zs|>;
    Export[FileNameJoin[{dir, "relaxed_eps" <> ToString[NumberForm[eps, {3, 2}]] <> ".csv"}], ReIm /@ zs];
    say["eps ", eps, ": map check max||w|-1| on the wall ", ScientificForm[chk, 2], "; flux check g0 = ", Chop[d["g0"], 10^-8],
      "; tip/flank density ratio ", NumberForm[tipRatio[d, zs], 4], " (uniform 1, sqrt(eps) ", NumberForm[Sqrt[eps], 4], ", conductor ", eps, "); depth ",
      NumberForm[Min[depths], 3], "-", NumberForm[Max[depths], 3], " cells"]],
  {eps, {1.0001, 1.2, 1.5, 2.0}}];
say["disk reference (analytic, N = 12, R = ", Sqrt[area/Pi], "): depth ", Sqrt[area/Pi] (1 - (11/35)^(1/24)), " cells"];

(* energy of a bulk test charge qt from the relaxed layer (+ images + boundary winding) and its own image, inside the 0.8 level *)
layerE[d_, zs_, qs_, z_, qt_] := Module[{w = d["w"][z], ws = d["w"] /@ zs},
  -2 qt Sum[qs[[i]] (Log[Abs[w - ws[[i]]]] + Log[Abs[1 - w Conjugate[ws[[i]]]]]), {i, Length[zs]}] - 2 Total[qs] qt d["H"][w]];
testE[d_, zs_, qs_, z_, qt_] := Module[{w = d["w"][z], ws = d["w"] /@ zs},
  -2 qt Sum[qs[[i]] (Log[Abs[w - ws[[i]]]] + Log[Abs[1 - w Conjugate[ws[[i]]]]]), {i, Length[zs]}]
  - qt^2 (Log[Abs[d["wp"][z]]] + Log[1 - Abs[w]^2]) - 2 (Total[qs] + qt) qt d["H"][w]];
Do[
  Module[{d = res[eps]["d"], zs = res[eps]["zs"], A, B, pts, vals, vp, vm, fig},
    A = d["A"]; B = d["B"];
    pts = Select[Flatten[Table[x + I y, {x, -A, A, A/18}, {y, -B, B, B/12}]], (Re[#]/A)^2 + (Im[#]/B)^2 < 0.8^2 &];
    vp = testE[d, zs, ConstantArray[1/2, nn], #, 1/2] & /@ pts;
    vm = testE[d, zs, ConstantArray[1/2, nn], #, -1/2] & /@ pts;
    Module[{lp = layerE[d, zs, ConstantArray[1/2, nn], #, 1/2] & /@ pts, inner},
      inner = (Re[#]/A)^2 + (Im[#]/B)^2 < 0.6^2 & /@ pts;
      say["eps ", NumberForm[eps, 4], ": LAYER part for a bulk +1/2 (pi K): range inside 0.8 level ", NumberForm[Max[lp] - Min[lp], 3],
        ", inside 0.6 level ", NumberForm[Max[Pick[lp, inner]] - Min[Pick[lp, inner]], 3], "; centre minus tip region ",
        NumberForm[Mean[Pick[lp, (Abs[Re[#]] < 0.2 A) & /@ pts]] - Mean[Pick[lp, (Abs[Re[#]] > A/Sqrt[2]) & /@ pts]], 3]]];
    say["eps ", NumberForm[eps, 4], ": energy range of a bulk test defect inside the 0.8 level (pi K): +1/2 ", NumberForm[Max[vp] - Min[vp], 3],
      ", -1/2 ", NumberForm[Max[vm] - Min[vm], 3], "; centre minus tip region (|x| > A/Sqrt[2]) for +1/2: ",
      NumberForm[Mean[Pick[vp, (Abs[Re[#]] < 0.2 A) & /@ pts]] - Mean[Pick[vp, (Abs[Re[#]] > A/Sqrt[2]) & /@ pts]], 3]];
    fig = Show[ListDensityPlot[Transpose[{Re[pts], Im[pts], vp - Mean[vp]}], ColorFunction -> "TemperatureMap", PlotRange -> {-0.5, 0.5},
        ClippingStyle -> Automatic, AspectRatio -> B/A, PlotLegends -> Automatic, FrameLabel -> {"x", "y"},
        PlotLabel -> Style["\[Epsilon] = " <> ToString[NumberForm[eps, 3]] <> ": energy of a bulk +1/2 (\[Pi]K, mean removed)", 10]],
      Graphics[{Black, PointSize[0.015], Point[ReIm /@ zs], Thin, Circle[{0, 0}, {A, B}]}], ImageSize -> 360];
    Export[FileNameJoin[{dir, "interior_potential_eps" <> ToString[NumberForm[eps, {3, 2}]] <> ".png"}], fig, ImageResolution -> 130]],
  {eps, {1.0001, 1.2, 1.5, 2.0}}];
Close[out];
