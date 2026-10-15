(* Ellipse shielding review (2026-10-15): maps to the 0.98 level; columns: layer, x, y, V(layer+probe)-V(layer) with flux = charges present (as before),
   probe self term, and the same V with the H (uniform-winding) coefficient set to the imposed winding K0 instead of the charge present. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "probe_maps_lib.wl"}]];
enQ[ws_, wps_, hs_, qs_, Q_] := Module[{n = Length[ws]},
  -Sum[2 qs[[i]] qs[[j]] (Log[Abs[ws[[i]] - ws[[j]]]] + Log[Abs[1 - ws[[i]] Conjugate[ws[[j]]]]]), {i, 1, n - 1}, {j, i + 1, n}]
  - Sum[qs[[i]]^2 (Log[Abs[wps[[i]]]] + Log[1 - Abs[ws[[i]]]^2]), {i, n}] - 2 Q Sum[qs[[i]] hs[[i]], {i, n}]];
dE = setupE[J["A2"], J["B2"]]; A = J["A2"]; B = J["B2"];
run2[layers_, name_, K0_] := Module[{grid, gw, gwp, gh, rows = {}},
  grid = Select[Flatten[Table[x + I y, {x, -A, A, 1.}, {y, -B, B, 1.}]], (Re[#]/A)^2 + (Im[#]/B)^2 < 0.98^2 &];
  gw = dE["w"] /@ grid; gwp = dE["wp"] /@ grid; gh = dE["H"] /@ gw;
  Do[Module[{L = layers[[l]], zs, qs, lw, lwp, lh, e0, e0k},
     zs = (#[[1]] + I #[[2]]) & /@ L; qs = #[[3]] & /@ L; lw = dE["w"] /@ zs; lwp = dE["wp"] /@ zs; lh = dE["H"] /@ lw;
     e0 = enQ[lw, lwp, lh, qs, Total[qs]]; e0k = enQ[lw, lwp, lh, qs, K0];
     Do[AppendTo[rows, {l, Re[grid[[p]]], Im[grid[[p]]],
         enQ[Append[lw, gw[[p]]], Append[lwp, gwp[[p]]], Append[lh, gh[[p]]], Append[qs, 1/2], Total[qs] + 1/2] - e0,
         enQ[{gw[[p]]}, {gwp[[p]]}, {gh[[p]]}, {1/2}, 1/2],
         enQ[Append[lw, gw[[p]]], Append[lwp, gwp[[p]]], Append[lh, gh[[p]]], Append[qs, 1/2], K0] - e0k}], {p, Length[grid]}]], {l, Length[layers]}];
  Export[FileNameJoin[{dir, "review_" <> name <> ".csv"}], rows]; Print[name]];
run2[J["ellipse_ideal"], "ellipse_ideal", 6]; run2[J["ellipse_meas"], "ellipse_meas", 6];
run2[J["mellipse_ideal"], "mellipse_ideal", -6]; run2[J["mellipse_meas"], "mellipse_meas", -6];
