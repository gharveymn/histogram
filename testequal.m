function idx = testequal(t1,t2)
  idx = 0;
  for i = 1:numel(t1)
    s1 = t1{i};
    s2 = t2{i};
    if(s1.counts != s2.counts || s1.edges != s2.edges || s1.binidx != s2.binidx)
      idx = i;
      return;
    endif
  endfor
endfunction
