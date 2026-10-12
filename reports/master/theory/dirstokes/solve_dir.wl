(* Stokes flow (eta = 1, no-slip ellipse/disk) of the measured active force density f = -div Q (prep_dir.py), per flow interval,
   separately for the wall band (within 4 cells of the wall) and the bulk. Force on the cell grid -> bilinear interpolation (zero outside).
   Velocities evaluated at the comparison cells. Output <tag>_uw.csv and <tag>_ub.csv: one row per interval, (ux, uy) per cell flattened.
   Run: wolframscript -file solve_dir.wl <tag> *)
Needs["NDSolve`FEM`"];
SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 2}];
dir = DirectoryName[$InputFileName]; tag = Last[$ScriptCommandLine];
meta = Import[FileNameJoin[{dir, tag <> "_meta.json"}], "RawJSON"];
{nint, nx, ny, X0, Y0, A, B, nc} = Lookup[meta, {"nint", "nx", "ny", "X0", "Y0", "A", "B", "ncells"}];
raw = BinaryReadList[FileNameJoin[{dir, tag <> "_force.bin"}], "Real64", ByteOrdering -> -1];
F = ArrayReshape[raw, {2, nint, nx, ny, 2}];
cells = Partition[BinaryReadList[FileNameJoin[{dir, tag <> "_cells.bin"}], "Real64", ByteOrdering -> -1], 2];
mesh = ToElementMesh[Disk[{0, 0}, {A, B}], MaxCellMeasure -> 0.35];
solve[g_] := Module[{fx, fy, sol},
  fx = ListInterpolation[g[[All, All, 1]], {{X0, X0 + nx - 1}, {Y0, Y0 + ny - 1}}, InterpolationOrder -> 1];
  fy = ListInterpolation[g[[All, All, 2]], {{X0, X0 + nx - 1}, {Y0, Y0 + ny - 1}}, InterpolationOrder -> 1];
  sol = Quiet[NDSolveValue[{{-Laplacian[u[x, y], {x, y}] + D[p[x, y], x] - fx[x, y], -Laplacian[v[x, y], {x, y}] + D[p[x, y], y] - fy[x, y],
       D[u[x, y], x] + D[v[x, y], y]} == {0, 0, 0}, DirichletCondition[{u[x, y] == 0, v[x, y] == 0}, True], DirichletCondition[p[x, y] == 0, x == A && y == 0]},
     {u, v}, Element[{x, y}, mesh], Method -> {"FiniteElement", "InterpolationOrder" -> {u -> 2, v -> 2, p -> 1}}]];
  Flatten[{sol[[1]] @@ #, sol[[2]] @@ #} & /@ cells]];
res = {{}, {}};
Do[Do[AppendTo[res[[r]], solve[F[[r, n]]]], {r, 2}]; If[Mod[n, 5] == 0, Print[tag, ": ", n, "/", nint]], {n, nint}];
Export[FileNameJoin[{dir, tag <> "_uw.csv"}], res[[1]]]; Export[FileNameJoin[{dir, tag <> "_ub.csv"}], res[[2]]];
Print["done ", tag];
