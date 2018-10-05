ts = histcountstestgen;
parfor i = 1:numel(ts)
  r1(i) = histcountstest(ts(i)); 
end
save('test.mat', '-v7');
