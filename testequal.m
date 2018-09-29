function idx = testequal(r1,r2)
  idx = 0;
  for i = 1:numel(r1)
    s1 = r1(i);
    s2 = r2(i);
    if(any(s1.c != s2.c) || any(s1.e != s2.e) || any(s1.b != s2.b))
      idx = i;
      return;
    endif
  endfor
endfunction
