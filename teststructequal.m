function retval = teststructequal(r1,r2)
  for i = 1:size(r1,1)
    rs1 = r1(i,:);
    rs2 = r2(i,:);
    for j = 1:numel(rs1)
      s1 = rs1(j);
      s2 = rs2(j);
      if(any(size(s1.c) != size(s2.c)) || any(size(s1.e) != size(s2.e)) || any(size(s1.b) != size(s2.b))...
         || any (abs (s1.c (! isnan (s1.c)) - s2.c (! isnan (s2.c))) > 100000000*eps) ...
         || any (abs (s1.e (! isnan (s1.e)) - s2.e (! isnan (s2.e))) > 100000000*eps) ...
         || any (abs (s1.b (! isnan (s1.b)) - s2.b (! isnan (s2.b))) > 100000000*eps) ...
         || ! strcmp(class(s1.c), class(s2.c))
         || ! strcmp(class(s1.e), class(s2.e))
         || ! strcmp(class(s1.b), class(s2.b)))
        retval = i;
        return;
      endif
    endfor
  endfor
  retval = 0;
endfunction
