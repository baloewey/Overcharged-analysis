(* Wall-layer tilt patterns in ellipses (area 2500): (1) static flows of a coherent and a split tilt; (2) feedback modes.
   Wall defects: N = 10 +1/2 evenly spaced in arclength at depth 1.5 (smeared point forces, sigma = 1.5, no-slip, eta = 1).
   (1) Force = cos(63 deg) n + s_j sin(63 deg) t (n outward normal, t ccw tangent); coherent: s_j = +1; split: s_j = +1 for x < 0, -1 for x > 0.
       Reported: rotation about the centroid of each half and of the whole cavity (f/eta).
   (2) M_ij = tangential (ccw) velocity at defect i from a unit tangential force at defect j (one solve per j). If the tilt follows the local
       tangential flow and the tangential force follows the tilt, M is the linear feedback operator; its eigenvectors are the self-reinforcing
       tilt patterns, the eigenvalues their gains. Mode classes: uniform (one vortex), dipolar (left/right opposite, two vortices), others.
   Run: wolframscript -file stokes_tilt_modes.wl *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{dir, "stokes_tilt_modes_results.txt"}]]; say[x__] := (WriteString[out, StringJoin[ToString /@ {x}], "\n"]; Print[x]);
sig = 1.5; nn = 10; depth = 1.5; tilt = 63 Degree; area = 2500.;
solve[A_, B_, forces_] := Module[{mesh = ToElementMesh[Disk[{0, 0}, {A, B}], MaxCellMeasure -> 0.35], fx, fy},
  fx = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[2, 1]] Exp[-((x - #[[1, 1]])^2 + (y - #[[1, 2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  fy = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[2, 2]] Exp[-((x - #[[1, 1]])^2 + (y - #[[1, 2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y], -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y],
      D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0}, DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == A && y == 0]},
   {u, v, p}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}]];
sites[A_, B_] := Module[{tt = Subdivide[0., 2 Pi, 4000], pts, s, P, ts},
  pts = {A Cos[#], B Sin[#]} & /@ tt; s = Prepend[Accumulate[Norm /@ Differences[pts]], 0.]; P = Last[s];
  ts = Interpolation[Transpose[{s, tt}], InterpolationOrder -> 1] /@ ((Range[nn] - 0.5) P/nn);
  Table[With[{t = ts[[j]]}, Module[{n = Normalize[{Cos[t]/A, Sin[t]/B}]}, {{A Cos[t], B Sin[t]} - depth n, n, {-n[[2]], n[[1]]}}]], {j, nn}]];
rotations[sol_, A_, B_] := Module[{pts = Select[Flatten[Table[{x, y}, {x, -A + 0.5, A - 0.5, 1.}, {y, -B + 0.5, B - 0.5, 1.}], 1], (#[[1]]/A)^2 + (#[[2]]/B)^2 < 0.95 &], rot},
  rot[ps_] := Module[{c = Mean[ps]}, Mean[(((#[[1]] - c[[1]]) sol[[2]] @@ # - (#[[2]] - c[[2]]) sol[[1]] @@ #)/Max[Norm[# - c], 0.5]) & /@ ps]];
  {rot[Select[pts, #[[1]] < 0 &]], rot[Select[pts, #[[1]] > 0 &]], rot[pts]}];
classify[vec_, st_] := Module[{sx = Sign[st[[All, 1, 1]]], u = Normalize[vec]},
  Which[Abs[Total[u]]/Sqrt[nn] > 0.8, "uniform", Abs[u . (-sx)]/Sqrt[nn] > 0.8, "dipolar (left/right)", True, "other"]];
Do[Module[{A = Sqrt[area eps/Pi], B, st, coh, spl, rc, rs, M, es, vs},
   B = A/eps; st = sites[A, B];
   If[MemberQ[{1.0001, 1.2, 1.5, 2.0}, eps],
    coh = solve[A, B, {#[[1]], Cos[tilt] #[[2]] + Sin[tilt] #[[3]]} & /@ st];
    spl = solve[A, B, {#[[1]], Cos[tilt] #[[2]] + If[#[[1, 1]] < 0, 1, -1] Sin[tilt] #[[3]]} & /@ st];
    rc = rotations[coh, A, B]; rs = rotations[spl, A, B];
    say["eps ", eps, ": coherent tilt -> rotation left/right/whole ", NumberForm[rc, 3], ";  split tilt -> ", NumberForm[rs, 3]]];
   M = Table[0., {nn}, {nn}];
   Do[Module[{sol = solve[A, B, {{st[[j, 1]], st[[j, 3]]}}]},
      Do[M[[i, j]] = {sol[[1]] @@ st[[i, 1]], sol[[2]] @@ st[[i, 1]]} . st[[i, 3]], {i, nn}]], {j, nn}];
   M = (M + Transpose[M])/2;
   {es, vs} = Eigensystem[M];
   say["eps ", eps, ": feedback gains (top 4) ", NumberForm[Take[es, 4], 3], "; modes: ", Table[classify[vs[[m]], st], {m, 4}],
     "; uniform gain ", NumberForm[Normalize[ConstantArray[1., nn]] . M . Normalize[ConstantArray[1., nn]], 3],
     ", dipolar gain ", NumberForm[With[{d = Normalize[-Sign[st[[All, 1, 1]]]]}, d . M . d], 3]]],
  {eps, {1.0001, 1.2, 1.35, 1.5, 1.65, 1.8, 2.0, 2.5}}];
Close[out];
