(* Coherently tilted wall layer in ellipses (area 2500): does it drive one vortex or two? 10 wall +1/2 evenly spaced in arclength at depth 1.5,
   head tilted +63 deg from the local outward normal (ccw sense everywhere). Halves x < 0 and x > 0: mean azimuthal velocity about each half's
   centroid (same sign: one vortex; opposite: two counter-rotating). Whole-cavity rotation: mean v_theta about the centre. Units f/eta. *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{dir, "stokes_ellipse_static_results.txt"}]]; say[x__] := (WriteString[out, StringJoin[ToString /@ {x}], "\n"]; Print[x]);
sig = 1.5;
solve[A_, B_, forces_] := Module[{mesh = ToElementMesh[Disk[{0, 0}, {A, B}], MaxCellMeasure -> 0.35], fx, fy},
  fx = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[2, 1]] Exp[-((x - #[[1, 1]])^2 + (y - #[[1, 2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  fy = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[2, 2]] Exp[-((x - #[[1, 1]])^2 + (y - #[[1, 2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y], -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y],
      D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0}, DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == A && y == 0]},
   {u, v, p}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}]];
layer[A_, B_, nn_, depth_, tilt_] := Module[{tt = Subdivide[0., 2 Pi, 4000], pts, s, P, targets, ts},
  pts = {A Cos[#], B Sin[#]} & /@ tt; s = Prepend[Accumulate[Norm /@ Differences[pts]], 0.]; P = Last[s];
  targets = (Range[nn] - 0.5) P/nn; ts = Interpolation[Transpose[{s, tt}], InterpolationOrder -> 1] /@ targets;
  Table[With[{t = ts[[j]]}, Module[{w = {A Cos[t], B Sin[t]}, nrm = Normalize[{Cos[t]/A, Sin[t]/B}], ang},
      ang = ArcTan @@ nrm + tilt; {w - depth nrm, {Cos[ang], Sin[ang]}}]], {j, nn}]];
halves[sol_, A_, B_] := Module[{pts = Select[Flatten[Table[{x, y}, {x, -A + 0.5, A - 0.5, 1.}, {y, -B + 0.5, B - 0.5, 1.}], 1], (#[[1]]/A)^2 + (#[[2]]/B)^2 < 0.95 &], L, Rr, rot},
  rot[ps_] := Module[{c = Mean[ps]}, Mean[(((#[[1]] - c[[1]]) sol[[2]] @@ # - (#[[2]] - c[[2]]) sol[[1]] @@ #)/Max[Norm[# - c], 0.5]) & /@ ps]];
  L = Select[pts, #[[1]] < 0 &]; Rr = Select[pts, #[[1]] > 0 &];
  {rot[L], rot[Rr], Mean[((#[[1]] sol[[2]] @@ # - #[[2]] sol[[1]] @@ #)/Max[Norm[#], 0.5]) & /@ pts]}];
figs = {};
Do[Module[{A = Sqrt[2500 eps/Pi], B, f, sol, h},
   B = A/eps; f = layer[A, B, 10, 1.5, 63 Degree]; sol = solve[A, B, f]; h = halves[sol, A, B];
   say["eps ", eps, ": rotation about left/right half centroids ", NumberForm[h[[1]], 3], " / ", NumberForm[h[[2]], 3], "; whole cavity ", NumberForm[h[[3]], 3], " (f/eta)"];
   AppendTo[figs, Show[StreamDensityPlot[{sol[[1]][x, y], sol[[2]][x, y]}, Element[{x, y}, Disk[{0, 0}, {A, B}]], ColorFunction -> "Rainbow",
       StreamStyle -> White, PlotLabel -> "\[Epsilon] = " <> ToString[eps], ImageSize -> 260],
      Graphics[{Black, Arrowheads[0.03], Arrow[{#[[1]], #[[1]] + 4 #[[2]]}] & /@ f}]]]],
  {eps, {1.0001, 1.2, 1.5, 2.0}}];
Export[FileNameJoin[{dir, "stokes_ellipse_static.png"}], GraphicsRow[figs, ImageSize -> 1100], ImageResolution -> 130];
Close[out];
