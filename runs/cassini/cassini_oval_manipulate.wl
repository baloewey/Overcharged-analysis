(* ::Package:: *)
(* Cassini-oval cavity explorer for H8 (figure-8 cavity).
   Oval d1*d2 = b^2, foci at (+-a, 0): ((x-a)^2+y^2)((x+a)^2+y^2) == b^4   (NIAMH-MPCD: "shape": 1, "R": b, "cassiniA": a).
   Controls: b/a (neck), half-length (sets the size: a = H/Sqrt[1+(b/a)^2]), winding k.
   b/a >= Sqrt[2]: convex; 1 < b/a < Sqrt[2]: concave neck; b/a = 1: lemniscate; b/a < 1: two separate ovals.
   Usage: Get["cassini_oval_manipulate.wl"] in a notebook (or paste), then evaluate  cassiniExplorer
*)

lActive = 22.2;   (* active length at zeta = 0.07, 200^2 box (RESEARCH.md) *)

cassiniPanel[ratio_?NumericQ, h_?NumericQ, k_Integer, showFoci_, showRef_, showDef_] :=
 Module[{a, b, ymax, curve, th, r2, pts, seg, cum, per, area, sOfTh, thOfS, defPts, waist, info,
         theta0, lobeArea, concave},
  a = h/Sqrt[1 + ratio^2]; b = ratio a;
  ymax = If[b >= Sqrt[2] a, Sqrt[b^2 - a^2], b^2/(2 a)];
  curve = ContourPlot[((x - a)^2 + y^2) ((x + a)^2 + y^2) == b^4,
    {x, -1.08 h, 1.08 h}, {y, -1.15 ymax - 2, 1.15 ymax + 2},
    PlotPoints -> 50, MaxRecursion -> 3, ContourStyle -> Directive[Thickness[0.006], GrayLevel[0.15]],
    Frame -> False, AspectRatio -> Automatic];
  If[b > a,
   (* single closed wall: polar form, arclength, area *)
   th = Subdivide[0., 2. Pi, 4000];
   r2 = a^2 Cos[2 th] + Sqrt[b^4 - a^4 Sin[2 th]^2];
   pts = Transpose[{Sqrt[r2] Cos[th], Sqrt[r2] Sin[th]}];
   seg = Norm /@ Differences[pts];
   cum = Prepend[Accumulate[seg], 0.];
   per = Last[cum];
   area = 0.5 Total[Most[r2]] (2. Pi/4000);
   thOfS = Interpolation[Transpose[{cum, th}], InterpolationOrder -> 1];
   defPts = Table[With[{t = thOfS[per (j + 0.5)/(2 Abs[k])]},
      Sqrt[a^2 Cos[2 t] + Sqrt[b^4 - a^4 Sin[2 t]^2]] {Cos[t], Sin[t]}], {j, 0, 2 Abs[k] - 1}];
   waist = 2 Sqrt[b^2 - a^2],
   (* b <= a: two ovals (or the lemniscate) *)
   theta0 = 0.5 ArcSin[Min[1., b^2/a^2]];
   lobeArea = NIntegrate[Sqrt[Max[0., b^4 - a^4 Sin[2 t]^2]], {t, -theta0, theta0}];
   area = 2 lobeArea; per = Missing["two walls"]; defPts = {}; waist = 0.
   ];
  concave = 1 < ratio < Sqrt[2];
  info = Grid[{
     {"a (cassiniA)", NumberForm[a, {5, 2}], "b (R)", NumberForm[b, {5, 2}]},
     {"half-length", Row[{NumberForm[h, {4, 1}], "  (", NumberForm[h/lActive, {3, 2}], " \[ScriptL]a)"}],
      "waist (full width)", If[b > a, Row[{NumberForm[waist, {4, 1}], "  (", NumberForm[waist/lActive, {3, 2}], " \[ScriptL]a)"}], "closed"]},
     {"area", NumberForm[area, {6, 0}], "particles (dens 20)", NumberForm[Round[20 area], DigitBlock -> 3]},
     {"perimeter", If[NumericQ[per], NumberForm[per, {5, 1}], "two walls"],
      "imposed density 2k/P", If[NumericQ[per], NumberForm[2. Abs[k]/per, {4, 3}], "\[Dash]"]},
     {"shape", Style[Which[ratio >= Sqrt[2], "convex", ratio > 1, "concave neck", ratio == 1, "lemniscate (neck closes)", True, "two separate ovals (topology change)"],
       Bold, If[concave, Darker[Orange], Black]], SpanFromLeft, SpanFromLeft}
     }, Alignment -> Left, Spacings -> {2, 0.6}];
  Column[{
    Show[curve,
     Graphics[{
       If[showRef, {Dashed, GrayLevel[0.55], Circle[{h - 1.434 lActive, 0}, 1.434 lActive], Circle[{-h + 1.434 lActive, 0}, 1.434 lActive]}, {}],
       If[showFoci, {RGBColor["#555a62"], PointSize[0.012], Point[{{-a, 0}, {a, 0}}]}, {}],
       If[showDef && defPts =!= {}, {RGBColor["#2a78d6"], PointSize[0.018], Point[defPts]}, {}],
       If[b > a, {RGBColor["#eb6834"], Thickness[0.004], Line[{{0, -Sqrt[b^2 - a^2]}, {0, Sqrt[b^2 - a^2]}}]}, {}]
       }],
     PlotRange -> {{-1.08 h, 1.08 h}, {-1.15 ymax - 2, 1.15 ymax + 2}}, ImageSize -> 560],
    info}, Spacings -> 1]
  ]

cassiniExplorer = Manipulate[
  cassiniPanel[ratio, h, k, showFoci, showRef, showDef],
  {{ratio, 1.3, "b/a (neck)"}, 0.9, 2.5, Appearance -> "Labeled"},
  {{h, 48, "half-length (cells)"}, 20, 80, 1, Appearance -> "Labeled"},
  {{k, 6, "winding k"}, 1, 12, 1, Appearance -> "Labeled"},
  {{showDef, True, "2k imposed +1/2 (uniform in arclength)"}, {True, False}},
  {{showFoci, True, "foci"}, {True, False}},
  {{showRef, False, "circles R = 1.43 \[ScriptL]a at the lobe ends"}, {True, False}},
  Delimiter,
  Style["Pilot P1 shapes: b/a = 1.6, 1.3, 1.1 at half-length 48 and 64; concave below b/a = \[Sqrt]2 \[TildeTilde] 1.41", 10],
  ControlPlacement -> Top, SaveDefinitions -> True];

cassiniExplorer

(* ---- Lobe coordinates: w = z^2 - a^2.  |w| = d1 d2 = c^2 (Cassini levels, c = 0: foci, c = a: lemniscate,
   c = b: wall), arg w = phi (orthogonal rectangular hyperbolas through the foci).
   Curve from a focus to the wall:  z(c; phi) = +-Sqrt[a^2 + c^2 E^(I phi)],  c in [0, b]  (+ right lobe, - left lobe). ---- *)

lobeCoordsPanel[ratio_?NumericQ, h_?NumericQ, nLevels_Integer, nRays_Integer] :=
 Module[{a, b, ymax, levels, rays, focusToWall},
  a = h/Sqrt[1 + ratio^2]; b = ratio a;
  ymax = If[b >= Sqrt[2] a, Sqrt[b^2 - a^2], b^2/(2 a)];
  levels = ContourPlot[Sqrt[Sqrt[(x - a)^2 + y^2] Sqrt[(x + a)^2 + y^2]], {x, -1.08 h, 1.08 h}, {y, -1.15 ymax - 2, 1.15 ymax + 2},
    Contours -> Table[b j/nLevels, {j, 1, nLevels - 1}], ContourShading -> None, PlotPoints -> 60,
    ContourStyle -> Directive[GrayLevel[0.65], Thin],
    RegionFunction -> Function[{x, y}, ((x - a)^2 + y^2) ((x + a)^2 + y^2) <= b^4], Frame -> False];
  focusToWall[s_, phi_] := ReIm[s Sqrt[a^2 + c^2 Exp[I phi]]];
  rays = ParametricPlot[Evaluate@Flatten[Table[{focusToWall[1, phi], focusToWall[-1, phi]}, {phi, 0., 2 Pi (1 - 1/nRays), 2 Pi/nRays}], 1],
    {c, 0, b}, PlotStyle -> Directive[RGBColor["#2a78d6"], Thickness[0.003]]];
  Show[levels, rays,
   ContourPlot[((x - a)^2 + y^2) ((x + a)^2 + y^2) == b^4, {x, -1.08 h, 1.08 h}, {y, -1.15 ymax - 2, 1.15 ymax + 2},
    ContourStyle -> Directive[Thickness[0.006], GrayLevel[0.15]], PlotPoints -> 50],
   ContourPlot[((x - a)^2 + y^2) ((x + a)^2 + y^2) == a^4, {x, -1.08 h, 1.08 h}, {y, -1.15 ymax - 2, 1.15 ymax + 2},
    ContourStyle -> Directive[RGBColor["#eb6834"], Dashed, Thickness[0.004]], PlotPoints -> 80, MaxRecursion -> 4],
   Graphics[{GrayLevel[0.2], PointSize[0.014], Point[{{-a, 0}, {a, 0}}]}],
   PlotRange -> {{-1.08 h, 1.08 h}, {-1.15 ymax - 2, 1.15 ymax + 2}}, AspectRatio -> Automatic, ImageSize -> 560,
   PlotLabel -> Style[Row[{"levels c = ", Row[{"b/", nLevels}], " steps (grey), c = a lemniscate (orange); curves focus \[Rule] wall, ", nRays, " angles (blue)"}], 10]]
  ]

lobeCoordsExplorer = Manipulate[
  lobeCoordsPanel[ratio, h, nLevels, nRays],
  {{ratio, 1.3, "b/a (neck)"}, 1.005, 2.5, Appearance -> "Labeled"},
  {{h, 48, "half-length (cells)"}, 20, 80, 1, Appearance -> "Labeled"},
  {{nLevels, 8, "Cassini levels"}, 2, 16, 1, Appearance -> "Labeled"},
  {{nRays, 16, "focus-to-wall curves per lobe"}, 4, 32, 1, Appearance -> "Labeled"},
  ControlPlacement -> Top, SaveDefinitions -> True];

lobeCoordsExplorer
