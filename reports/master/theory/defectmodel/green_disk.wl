(* No-slip Stokes Green's function of the unit disk (eta = 1), for the defect Langevin model. In 2D, u = (F/eta) g(z/R; z_s/R) is
   independent of the disk radius, so one table serves every R. Point force (Gaussian smearing sig = 0.012, ~0.4 cells at R = 32) at
   z_s = (r_s, 0), r_s = 0..0.96 (49 values), along x and along y. Velocity sampled on a polar grid r = 0..1 (101), phi = 0..2 pi (181).
   Output green_disk.bin (Real64, order [rs][dir][r][phi][component]) and green_disk_meta.json. Run: wolframscript -file green_disk.wl *)
Needs["NDSolve`FEM`"];
dir = DirectoryName[$InputFileName]; sig = 0.012;
mesh = ToElementMesh[Disk[{0, 0}, 1], MaxCellMeasure -> 0.0003, "MaxBoundaryCellMeasure" -> 0.01];
Print["mesh nodes ", Length[mesh["Coordinates"]]];
rsList = Range[0, 0.96, 0.02]; rr = Subdivide[0., 1., 100]; pp = Subdivide[0., 2 Pi, 180];
solve[xs_, d_] := Module[{fx, fy, sol},
  fx = Function[{x, y}, Evaluate[d[[1]] Exp[-((x - xs)^2 + y^2)/(2 sig^2)]/(2 Pi sig^2)]];
  fy = Function[{x, y}, Evaluate[d[[2]] Exp[-((x - xs)^2 + y^2)/(2 sig^2)]/(2 Pi sig^2)]];
  sol = Quiet[NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y], -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y],
       D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0}, DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == 1 && y == 0]},
     {u, v}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}]];
  Table[With[{x = Min[r, 0.9999] Cos[ph], y = Min[r, 0.9999] Sin[ph]}, {sol[[1]][x, y], sol[[2]][x, y]}], {r, rr}, {ph, pp}]];
out = Table[Print["r_s = ", rs]; Table[solve[rs, d], {d, {{1, 0}, {0, 1}}}], {rs, rsList}];
BinaryWrite[FileNameJoin[{dir, "green_disk.bin"}], Flatten[out], "Real64", ByteOrdering -> -1]; Close[FileNameJoin[{dir, "green_disk.bin"}]];
Export[FileNameJoin[{dir, "green_disk_meta.json"}], <|"rs" -> rsList, "nr" -> Length[rr], "nphi" -> Length[pp], "sig" -> sig|>, "RawJSON"];
Print["done"];
