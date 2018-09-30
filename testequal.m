function idx = testequal(r1,r2)
  idx = 0;
  for j = 1:numel(r1)
    currr1 = r1{j};
    currr2 = r2{j};
    for i = 1:numel(currr1)
      s1 = currr1(i);
      s2 = currr2(i);
      if(any(abs(s1.c-s2.c) > 100000000*eps) || any(abs(s1.e-s2.e) > 100000000*eps) || any(abs(s1.b-s2.b) > 100000000*eps)
         || any(size(s1.c) != size(s2.c)) || any(size(s1.e) != size(s2.e)) || any(size(s1.b) != size(s2.b)))
        idx = i;
        return;
      endif
    endfor
  endfor
endfunction
