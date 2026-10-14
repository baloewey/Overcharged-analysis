Get[FileNameJoin[{DirectoryName[$InputFileName], "probe_maps_lib.wl"}]];
dC = setupC[J["R1"]]; dE = setupE[J["A2"], J["B2"]];
maps[dC, J["mcircle_ideal"], {J["R1"], J["R1"]}, "mcircle_ideal"]; maps[dC, J["mcircle_meas"], {J["R1"], J["R1"]}, "mcircle_meas"];
maps[dE, J["mellipse_ideal"], {J["A2"], J["B2"]}, "mellipse_ideal"]; maps[dE, J["mellipse_meas"], {J["A2"], J["B2"]}, "mellipse_meas"];
