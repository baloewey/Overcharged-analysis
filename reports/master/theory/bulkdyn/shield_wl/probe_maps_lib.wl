(* Potential felt by one bulk +1/2 from the wall layer alone (layer charges + their wall images + the boundary winding), one-constant, units of pi K.
   Same construction as the master document (Sec. 'Elastic theory of the wall layer: shielding', theory/shielding_ellipse.wl) and the Japan poster:
   V(z) = E(layer + probe at z) - E(layer), with E from the uniform-flux Neumann Green's function (circle: images at R^2/conj(z); ellipse: conformal map).
   Input: layers.json (python: measured wall layers from the fine-output runs, and idealised full layers). Output: maps_<case>.csv (x, y, V per layer).
   Run: wolframscript -file probe_maps.wl *)
dir = DirectoryName[$InputFileName];
J = Import[FileNameJoin[{dir, "layers.json"}], "RawJSON"];

setupE[A_, B_] := Module[{c, q, m, K, wf, wp, zf, P, nF = 256, al, zs, g, cf},
  c = Sqrt[A^2 - B^2]; q = ((A - B)/(A + B))^2; m = InverseEllipticNomeQ[q]; K = EllipticK[m];
  wf = Function[z, m^(1/4) JacobiSN[(2 K/Pi) ArcSin[z/c], m]];
  wp = Function[z, Evaluate[D[m^(1/4) JacobiSN[(2 K/Pi) ArcSin[zz/c], m], zz] /. zz -> z]];
  zf = Function[w, c Sin[(Pi/(2 K)) InverseJacobiSN[w/m^(1/4), m]]];
  P = 4 A EllipticE[1 - B^2/A^2];
  al = Table[2 Pi j/nF, {j, 0, nF - 1}]; zs = zf /@ Exp[I al];
  g = (2 Pi/P)/Abs[wp /@ zs] - 1; cf = Fourier[g, FourierParameters -> {-1, -1}];
  <|"w" -> wf, "wp" -> wp, "H" -> Function[w, Re[Sum[2 cf[[n + 1]] w^n/n, {n, 1, nF/2 - 1}]]], "inside" -> Function[z, (Re[z]/A)^2 + (Im[z]/B)^2 < 0.9^2]|>];
setupC[R_] := <|"w" -> Function[z, z/R], "wp" -> Function[z, 1/R], "H" -> Function[w, 0.], "inside" -> Function[z, Abs[z] < 0.9 R]|>;

(* energy from precomputed w, w', H of every charge *)
en[ws_, wps_, hs_, qs_] := Module[{n = Length[ws], Q = Total[qs]},
  -Sum[2 qs[[i]] qs[[j]] (Log[Abs[ws[[i]] - ws[[j]]]] + Log[Abs[1 - ws[[i]] Conjugate[ws[[j]]]]]), {i, 1, n - 1}, {j, i + 1, n}]
  - Sum[qs[[i]]^2 (Log[Abs[wps[[i]]]] + Log[1 - Abs[ws[[i]]]^2]), {i, n}] - 2 Q Sum[qs[[i]] hs[[i]], {i, n}]];

maps[d_, layers_, ext_, name_] := Module[{grid, gw, gwp, gh, rows = {}},
  grid = Select[Flatten[Table[x + I y, {x, -ext[[1]], ext[[1]], 1.}, {y, -ext[[2]], ext[[2]], 1.}]], d["inside"]];
  gw = d["w"] /@ grid; gwp = d["wp"] /@ grid; gh = d["H"] /@ gw;
  Do[Module[{L = layers[[l]], zs, qs, lw, lwp, lh, e0},
     zs = (#[[1]] + I #[[2]]) & /@ L; qs = #[[3]] & /@ L; lw = d["w"] /@ zs; lwp = d["wp"] /@ zs; lh = d["H"] /@ lw; e0 = en[lw, lwp, lh, qs];
     Do[AppendTo[rows, {l, Re[grid[[p]]], Im[grid[[p]]],
         en[Append[lw, gw[[p]]], Append[lwp, gwp[[p]]], Append[lh, gh[[p]]], Append[qs, 1/2]] - e0,
         en[{gw[[p]]}, {gwp[[p]]}, {gh[[p]]}, {1/2}]}], {p, Length[grid]}]],   (* 5th column: the probe alone (its own wall images) *)
   {l, Length[layers]}];
  Export[FileNameJoin[{dir, "maps_" <> name <> ".csv"}], rows]; Print[name, ": ", Length[layers], " layers x ", Length[grid], " points"]];

