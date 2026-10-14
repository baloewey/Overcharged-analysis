(* Elastic hold on each wall +1/2 of a layer: energy cost to move it 6 cells inward along the ellipse normal, the rest of the layer fixed (layer alone,
   images + boundary winding, one-constant, pi K units). Output: binding_<case>.csv (layer, j, x, y, dE). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "probe_maps_lib.wl"}]];
dE = setupE[J["A2"], J["B2"]]; A = J["A2"]; B = J["B2"];
run[layers_, name_, nws_:None] := Module[{rows = {}},
  Do[Module[{L = layers[[l]], zs, qs, e0},
     zs = (#[[1]] + I #[[2]]) & /@ L; qs = #[[3]] & /@ L;
     e0 = en[dE["w"] /@ zs, dE["wp"] /@ zs, dE["H"] /@ (dE["w"] /@ zs), qs];
     Do[Module[{z = zs[[j]], n, z2, zz},  (* only the first nw entries are wall defects *)
        n = -(Re[z]/A^2 + I Im[z]/B^2); n = n/Abs[n]; z2 = z + 6 n; zz = ReplacePart[zs, j -> z2];
        AppendTo[rows, {l, j, Re[z], Im[z], en[dE["w"] /@ zz, dE["wp"] /@ zz, dE["H"] /@ (dE["w"] /@ zz), qs] - e0}]], {j, If[nws === None, Length[zs], nws[[l]]]}]], {l, Length[layers]}];
  Export[FileNameJoin[{dir, "binding_" <> name <> ".csv"}], rows]; Print[name]];
run[J["ellipse_full"], "ellipse_full", J["ellipse_full_nwall"]];
