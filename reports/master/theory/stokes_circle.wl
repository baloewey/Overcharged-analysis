(* Static Stokes flows of defect forces in a no-slip disk (B. Loewe's idea, 2026-10-06).
   -eta Laplacian u + grad p = sum of smeared point forces, div u = 0, u = 0 on the wall. eta = 1, unit force per defect,
   so velocities are in units of f/eta. +1/2: one point force along its orientation (head). -1/2: three point forces at 120 deg,
   pointing inwards, each at distance a from the centre (Norton et al. 2018), magnitude f each. Smearing: Gaussian, width sig.
   Circulation as measured: spatial mean of v_theta over r < R - 1. Profiles: azimuthal average of v_theta in 2-cell bins.
   Run: wolframscript -file stokes_circle.wl  -> stokes_circle_results.txt, stokes_circle_profiles.csv, stokes_circle_fields.png *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{dir, "stokes_circle_results.txt"}]];
say[x__] := (WriteString[out, StringJoin[ToString /@ {x}], "\n"]; Print[x]);

R = 31.83; sig = 1.5; a3 = 1.5; tilt = 63 Degree;
mesh = ToElementMesh[Disk[{0, 0}, R], MaxCellMeasure -> 0.35];
say["mesh: ", Length[mesh["Coordinates"]], " nodes"];

(* force density from a list of {position {x,y}, vector {fx,fy}} *)
fdens[forces_] := Function[{x, y}, Evaluate[Total[(#[[2]] Exp[-((x - #[[1, 1]])^2 + (y - #[[1, 2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
solve[forces_] := Module[{fd = fdens[forces], fx, fy, sol},
  fx[x_, y_] := fd[x, y][[1]]; fy[x_, y_] := fd[x, y][[2]];
  sol = NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y],
       -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y], D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0},
     DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == R && y == 0]},
    {u, v, p}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}];
  sol];
vtheta[sol_, x_, y_] := (x sol[[2]][x, y] - y sol[[1]][x, y])/Sqrt[x^2 + y^2];
circulation[sol_] := Module[{pts = Select[Flatten[Table[{x, y}, {x, -R + 0.5, R - 0.5, 1.}, {y, -R + 0.5, R - 0.5, 1.}], 1], Norm[#] < R - 1 &]},
  Mean[vtheta[sol, #[[1]], #[[2]]] & /@ pts]];
profile[sol_] := Table[Mean[Table[vtheta[sol, r Cos[t], r Sin[t]], {t, 0, 2 Pi - 0.01, 2 Pi/72}]], {r, 1, R - 1, 2}];

(* wall layer: N +1/2 at the active depth (bulk charge 1 at the centre, Sec. theory) + 0.61-cell offset; head tilted by +63 deg from the outward normal (ccw) *)
wallLayer[nn_] := Module[{rw = R ((nn - 1 + 4)/(3 nn - 1 + 4))^(1/(2 nn)) - 0.61},
  Table[With[{ph = 2 Pi (j - 0.5)/nn}, {rw {Cos[ph], Sin[ph]}, {Cos[ph + tilt], Sin[ph + tilt]}}], {j, nn}]];
(* bulk +1: two co-rotating +1/2 at radius rb, heads along +theta (ccw) *)
bulkPair[rb_] := Table[With[{ph = Pi j}, {rb {Cos[ph], Sin[ph]}, {-Sin[ph], Cos[ph]}}], {j, 0, 1}];
(* -1/2 at position c with orientation angle o: three inward forces *)
tripole[c_, o_] := Table[With[{t = o + 2 Pi j/3}, {c + a3 {Cos[t], Sin[t]}, -{Cos[t], Sin[t]}}], {j, 0, 2}];

rows = {}; profs = {};
Do[Module[{s = solve[wallLayer[nn]], c},
   c = circulation[s]; AppendTo[rows, {"wall", nn, c}]; AppendTo[profs, Prepend[profile[s], nn]];
   say["wall layer N = ", nn, " (depth ", NumberForm[R - Norm[wallLayer[nn][[1, 1]]], 3], "): circulation ", NumberForm[c, 4], " f/eta, per wall defect ", NumberForm[c/nn, 4]]],
  {nn, {2, 6, 10, 14, 18}}];
Do[Module[{s = solve[bulkPair[rb R]], c},
   c = circulation[s]; AppendTo[rows, {"bulk", rb, c}]; AppendTo[profs, Prepend[profile[s], 100 + rb]];
   say["bulk pair at r = ", rb, " R: circulation ", NumberForm[c, 4], " f/eta"]],
  {rb, {0.3, 0.5, 0.6, 0.7, 0.8}}];
Module[{s = solve[tripole[{0, 0}, 0]]}, say["one -1/2 tripole at the centre: circulation ", NumberForm[circulation[s], 3],
   ", max |u| ", NumberForm[Max[Norm[{s[[1]][#[[1]], #[[2]]], s[[2]][#[[1]], #[[2]]]}] & /@ RandomPoint[Disk[{0, 0}, R - 1], 400]], 3], " f/eta"]];
(* single force at depth h from the wall, tangential: hydrodynamic screening by the no-slip wall *)
Do[Module[{s = solve[{{(R - h) {1, 0}, {0, 1}}}]}, say["single tangential force at depth ", h, ": circulation ", NumberForm[circulation[s], 4],
    ", velocity at the force ", NumberForm[s[[2]][R - h, 0], 4]]], {h, {1, 2, 4, 8, 16, 31.83}}];
Export[FileNameJoin[{dir, "stokes_circle_profiles.csv"}], profs];
(* fields: wall layer N = 10 and the bulk pair at 0.6 R *)
Module[{sw = solve[wallLayer[10]], sb = solve[bulkPair[0.6 R]], plt},
  plt[s_, ttl_, fr_] := Show[StreamDensityPlot[{s[[1]][x, y], s[[2]][x, y]}, Element[{x, y}, Disk[{0, 0}, R]],
      ColorFunction -> "Rainbow", StreamStyle -> White, PlotLegends -> Automatic, PlotLabel -> ttl, ImageSize -> 340],
    Graphics[{Black, Arrowheads[0.03], Arrow[{#[[1]], #[[1]] + 4 #[[2]]}] & /@ fr}]];
  Export[FileNameJoin[{dir, "stokes_circle_fields.png"}],
   GraphicsRow[{plt[sw, "10 wall +1/2, tilt 63\[Degree] (speed, f/\[Eta])", wallLayer[10]], plt[sb, "bulk pair at 0.6R (speed, f/\[Eta])", bulkPair[0.6 R]]}, ImageSize -> 760],
   ImageResolution -> 130]];
Close[out];
