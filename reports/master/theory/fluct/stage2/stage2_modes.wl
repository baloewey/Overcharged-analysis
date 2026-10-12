(* Stage 2: Stokes feedback matrix M_ij (tangential velocity at site i from a unit tangential force at site j, eta = 1) for each configuration
   in cfg_<tag>.json (make_configs.py). Same solver as stokes_tilt_modes.wl (smeared force sigma = 1.5, no-slip), mesh built once per file.
   Run: wolframscript -file stage2_modes.wl e15   -> M_e15.json (list of {id, M}) *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName]; tag = Last[$ScriptCommandLine];
cfgs = Import[FileNameJoin[{dir, "cfg_" <> tag <> ".json"}], "RawJSON"];
sig = 1.5; {A, B} = {cfgs[[1]]["A"], cfgs[[1]]["B"]};
meshes = <||>;
meshFor[a_, b_] := Lookup[meshes, Key[{a, b}], meshes[{a, b}] = ToElementMesh[Disk[{0, 0}, {a, b}], MaxCellMeasure -> 0.35]];
solve[mesh_, a_, pos_, vec_] := Module[{fx, fy},
  fx = Compile[{{x, _Real}, {y, _Real}}, Evaluate[vec[[1]] Exp[-((x - pos[[1]])^2 + (y - pos[[2]])^2)/(2 sig^2)]/(2 Pi sig^2)]];
  fy = Compile[{{x, _Real}, {y, _Real}}, Evaluate[vec[[2]] Exp[-((x - pos[[1]])^2 + (y - pos[[2]])^2)/(2 sig^2)]/(2 Pi sig^2)]];
  Quiet[NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y], -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y],
      D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0}, DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == a && y == 0]},
   {u, v}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}]]];
res = {};
Do[Module[{c = cfgs[[n]], mesh, st, nn, M},
   mesh = meshFor[c["A"], c["B"]]; st = c["sites"]; nn = Length[st];
   M = Table[0., {nn}, {nn}];
   Do[Module[{sol = solve[mesh, c["A"], st[[j, 1 ;; 2]], st[[j, 3 ;; 4]]]},
      Do[M[[i, j]] = {sol[[1]] @@ st[[i, 1 ;; 2]], sol[[2]] @@ st[[i, 1 ;; 2]]} . st[[i, 3 ;; 4]], {i, nn}]], {j, nn}];
   AppendTo[res, <|"id" -> c["id"], "M" -> M|>];
   If[Mod[n, 10] == 0, Print[tag, ": ", n, "/", Length[cfgs]];
     Export[FileNameJoin[{dir, "M_" <> tag <> ".json"}], res, "RawJSON"]]], {n, Length[cfgs]}];
Export[FileNameJoin[{dir, "M_" <> tag <> ".json"}], res, "RawJSON"];
Print["done ", tag];
