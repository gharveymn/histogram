load m6.mat
pkg load parallel
disp ("loaded");
r1 = pararrayfun (8, @histcountstest, ts);
disp ("received test results");
[x, fts, fr1, fr2] = testequal (ts, r, r1);
if (isempty (fts))
  disp ("test passed");
else
  disp ("test failed");
endif
