(* wrapper: limit threads per kernel (8 kernels run side by side), then run stage2_modes.wl on the tag given on the command line *)
SetSystemOptions["MKLThreads" -> 2, "ParallelOptions" -> {"ParallelThreadNumber" -> 2}];
Get[FileNameJoin[{DirectoryName[$InputFileName], "stage2_modes.wl"}]];
