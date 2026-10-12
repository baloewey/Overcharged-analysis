(* For each flow interval: Stokes flows (eta = 1, unit force per defect, smeared Gaussian width 1.5) of
   (1) wall +1/2 forces along the defect angle, (2) bulk +1/2 forces along the defect angle, (3) -1/2 tripoles (three inward forces at
   angle + 2 pi j/3, distance 1.5); no-slip disk. Velocities at the cell centres -> CSV per run and class.
   Run: wolframscript -file solve.wl k02   (argument: run tag) *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName]; tag = SelectFirst[$ScriptCommandLine, StringMatchQ[#, ("k" | "e") ~~ DigitCharacter ..] &];
dat = Import[FileNameJoin[{dir, tag <> ".json"}], "RawJSON"];
cells = dat["cells"]; sig = 1.5; a3 = 1.5;
If[KeyExistsQ[dat, "R"], {Ax, Bx} = {dat["R"], dat["R"]}, {Ax, Bx} = {dat["A"], dat["B"]}]; R = Ax;
mesh = ToElementMesh[Disk[{0, 0}, {Ax, Bx}], MaxCellMeasure -> 0.35];
solve[forces_] := Module[{fx, fy, sol},
  If[forces === {}, Return[ConstantArray[{0., 0.}, Length[cells]]]];
  fx = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[3]] #[[5]] Exp[-((x - #[[1]])^2 + (y - #[[2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  fy = Compile[{{x, _Real}, {y, _Real}}, Evaluate[Total[(#[[4]] #[[5]] Exp[-((x - #[[1]])^2 + (y - #[[2]])^2)/(2 sig^2)]/(2 Pi sig^2)) & /@ forces]]];
  sol = NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y],
       -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y], D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0},
     DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == R && y == 0]},
    {u, v, p}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}];
  {sol[[1]] @@ #, sol[[2]] @@ #} & /@ cells];
tripoleForces[ms_] := Flatten[Function[m, Table[With[{t = m[[3]] + 2 Pi j/3}, {m[[1]] + a3 Cos[t], m[[2]] + a3 Sin[t], -Cos[t], -Sin[t], m[[4]]}], {j, 0, 2}]] /@ ms, 1];
res = <|"wall" -> {}, "bulk" -> {}, "minus" -> {}|>;
onlyMinus = MemberQ[$ScriptCommandLine, "minus"];
Do[Module[{f = dat["frames"][[i]]},
   If[!onlyMinus, AppendTo[res["wall"], Flatten[solve[f["wallp"]]]]; AppendTo[res["bulk"], Flatten[solve[f["bulkp"]]]]];
   AppendTo[res["minus"], Flatten[solve[If[f["minus"] === {}, {}, tripoleForces[f["minus"]]]]]];
   If[Mod[i, 10] == 0, Print[tag, ": ", i, " frames"]]], {i, Length[dat["frames"]]}];
Do[Export[FileNameJoin[{dir, tag <> "_" <> c <> ".csv"}], res[c]], {c, If[onlyMinus, {"minus"}, {"wall", "bulk", "minus"}]}];
Print["done ", tag];
