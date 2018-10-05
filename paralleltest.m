clear
addpath tests
if (exist ("tests/test.mat", "file") == 0)
  matlabscript = ["cd('" pwd "/tests'); ts = histcountstestgen; " ...
                   "r1 = histcountstest(ts); save('test.mat', '-v7'); exit"];
  system(["matlab -nodesktop -r \"" matlabscript "\""]);
endif
load tests/test.mat
pkg load parallel
disp ("loaded");
r2 = pararrayfun (8, @histcountstest, ts);
disp ("received test results");
retvals = parcellfun(8, @teststructequal, r1, r2);
if (any (retvals))
  disp ("test failed");
else
  disp ("test passed");
endif
