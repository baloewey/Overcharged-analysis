(* Shielding in a disk: one-constant nematic, strong anchoring (fixed director angle on the wall = Neumann condition for the
   conjugate field psi, same-sign images at R^2/zbar). Energies in units of pi K.
   E = - Sum_{i != j} q_i q_j Log|z_i - z_j| - Sum_{i,j} q_i q_j Log|R^2 - z_i Conjugate[z_j]|        (core terms dropped)
   Run: wolframscript -file shielding_disk.wl   (writes PNGs and shielding_disk_results.txt next to this file) *)

dir = DirectoryName[$InputFileName];
out = OpenWrite[FileNameJoin[{dir, "shielding_disk_results.txt"}]];
say[x__] := (WriteString[out, StringJoin[ToString /@ {x}], "\n"]; Print[x]);

(* 1. ring identities *)
say["1. Sum_j Log|z - r w^j| = Log|z^N - r^N| (N = 7, random z):  ",
  N[Sum[Log[Abs[z - r Exp[2 Pi I j/7]]], {j, 0, 6}] - Log[Abs[z^7 - r^7]] /. {z -> 0.31 + 0.47 I, r -> 0.9}]];
say["   angular average of Log|1 - a Exp[I t]| (a = 0.3, 0.7, 0.95; exact value 0 for a < 1): ",
  Table[NIntegrate[Log[Abs[1 - a Exp[I t]]], {t, 0, 2 Pi}]/(2 Pi), {a, {0.3, 0.7, 0.95}}]];

(* 2. potential of a ring (N charges qw at radius rw, angles 2 pi j/N) and its images on a test charge qt at r Exp[I phi] *)
Vring[r_, phi_, nn_, rw_, R_, qt_, qw_] := -2 qt qw (nn Log[rw] + Log[Abs[1 - (r/rw)^nn Exp[I nn phi]]] +
    2 nn Log[R] + Log[Abs[1 - (r rw/R^2)^nn Exp[I nn phi]]]);
Vself[r_, R_, qt_] := -qt^2 Log[R^2 - r^2];
say["2. ring potential = const + corrugation of relative size (r/rw)^N; the angular average is independent of r and N, so the bulk",
    " does not feel the winding (exact for a circle)."];

(* 3. equilibrium radius of a passive ring of N charges q *)
Ering[r_, nn_, R_, q_] := -q^2 nn ((nn - 1) Log[r] + Log[R^(2 nn) - r^(2 nn)]);
rstar = R ((n - 1)/(3 n - 1))^(1/(2 n));
say["3. claimed minimum r* = R ((N-1)/(3N-1))^(1/2N): dE/dr at r* simplifies to ",
  FullSimplify[D[Ering[r, n, R, q], r] /. r -> rstar, n > 1 && R > 0 && q > 0]];
say["   second derivative at r* (N = 12, R = 1): ", N[D[Ering[r, 12, 1, 1/2], {r, 2}] /. r -> (11/35)^(1/24)], " (> 0: minimum)"];
say["   large-N depth: d/spacing -> ", N[Log[3]/(4 Pi)], "  (spacing = 2 pi R / N)"];
(* with a bulk charge Qb at the centre (active layer): extra term -2 q Qb N Log[r] *)
Eact[r_, nn_, R_, q_, Qb_] := Ering[r, nn, R, q] - 2 q Qb nn Log[r];
rstara = R ((n - 1 + 2 Qb/q)/(3 n - 1 + 2 Qb/q))^(1/(2 n));
say["   with a central bulk charge Qb, claimed r* = R ((N-1+2Qb/q)/(3N-1+2Qb/q))^(1/2N): dE/dr at r* = ",
  FullSimplify[D[Eact[r, n, R, q, Qb], r] /. r -> rstara, n > 1 && R > 0 && q > 0 && Qb >= 0]];

(* 4. approach of a bulk +1/2 to a wall -1/2 (negative winding) vs to a free -1/2 *)
(* circle R = 31.83 (1.43 l_a), k = -6: active wall lattice 2|k|+2 = 14 sites; ring depth from the measured passive law + offset *)
R0 = 31.83; nW = 14; spacing = 2 Pi R0/nW;
rw = R0 - (0.61 + R0 (1 - ((nW - 1)/(3 nW - 1))^(1/(2 nW))));
say["4. R = ", R0, ", N = ", nW, " wall -1/2, spacing ", N[spacing], ", ring at depth ", N[R0 - rw], " cells"];
(* test +1/2 on the radius through a wall defect (phi = 0), distance d from it *)
Vwall[d_] := Vring[rw - d, 0, nW, rw, R0, 1/2, -1/2] + Vself[rw - d, R0, 1/2] - (Vring[rw - d, 0, nW, rw, R0, 1/2, -1/2] /. {} );
Vline[d_] := Vring[rw - d, 0, nW, rw, R0, 1/2, -1/2] + Vself[rw - d, R0, 1/2];
Vmid[d_] := Vring[rw - d, Pi/nW, nW, rw, R0, 1/2, -1/2] + Vself[rw - d, R0, 1/2];     (* between two wall defects *)
Vfree[d_] := -2 (1/2) (-1/2) Log[d];                                                   (* free pair, attraction *)
W[f_, d1_, d2_] := N[f[d2] - f[d1]];  (* energy released moving from d2 to d1 (positive = attraction does work) *)
say["   energy released bringing the +1/2 from d = 4 to d = 1 cells (units of pi K):"];
say["     free -1/2 pair:                 ", W[Vfree, 1, 4]];
say["     towards a wall -1/2 (radial):   ", W[Vline, 1, 4], "   ratio to free: ", W[Vline, 1, 4]/W[Vfree, 1, 4]];
say["     from d = 8 to d = 4 (approach): free ", W[Vfree, 4, 8], ", wall ", W[Vline, 4, 8]];
say["   force at d = 4 (towards the defect > 0): free ", N[-D[Vfree[d], d] /. d -> 4], ", wall ", N[-D[Vline[d], d] /. d -> 4]];

pl1 = Plot[{Vfree[d] - Vfree[8], Vline[d] - Vline[8], Vmid[d] - Vmid[8]}, {d, 0.8, 12},
  PlotLegends -> {"free -1/2", "towards a wall -1/2", "between two wall -1/2"}, Frame -> True,
  FrameLabel -> {"distance d of the +1/2 from the wall layer (cells; radial path)", "energy (\[Pi]K), zero at d = 8"},
  PlotLabel -> Style["Bulk +1/2 approaching the wall -1/2 layer (k = -6, R = 1.43 \[ScriptL]a, 14 wall sites)", 11],
  PlotStyle -> {Directive[Gray, Dashed], Directive[RGBColor["#2a78d6"], Thick], Directive[RGBColor["#eb6834"], Thick]}, ImageSize -> 620];
Export[FileNameJoin[{dir, "approach_energy.png"}], pl1, ImageResolution -> 150];

(* 5. interior potential of the ring: map, showing the flat (shielded) interior *)
pl2 = DensityPlot[If[x^2 + y^2 < R0^2, Vring[Sqrt[x^2 + y^2], Arg[x + I y], nW, rw, R0, 1/2, -1/2] + Vself[Sqrt[x^2 + y^2], R0, 1/2], None],
  {x, -R0, R0}, {y, -R0, R0}, PlotPoints -> 160, ColorFunction -> "TemperatureMap", PlotRange -> {-2, 1}, ClippingStyle -> Automatic,
  PlotLegends -> Automatic, FrameLabel -> {"x", "y"}, PlotLabel -> Style["Energy of a bulk +1/2 (\[Pi]K): 14 wall -1/2 + images + its own image", 10],
  ImageSize -> 420];
Export[FileNameJoin[{dir, "interior_potential_disk.png"}], pl2, ImageResolution -> 150];
say["5. interior potential: variation over r < 0.7 R: ",
  N[(Max[#] - Min[#]) &@Flatten@Table[Vring[r, p, nW, rw, R0, 1/2, -1/2], {r, 0, 0.7 R0, 0.7 R0/20}, {p, 0, 2 Pi, 2 Pi/40}]], " (pi K)"];
Close[out];
