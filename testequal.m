function [x, fts, fr1, fr2] = testequal(ts,r1,r2)
  x   = [];
  fts = [];
  fr1 = [];
  fr2 = [];
  for j = 1:numel(r1)
    for i = 1:numel(r1{j})
      s1 = r1{j}(i);
      s2 = r2{j}(i);
      if(any(size(s1.c) != size(s2.c)) || any(size(s1.e) != size(s2.e)) || any(size(s1.b) != size(s2.b))...
         || any (abs (s1.c (! isnan (s1.c)) - s2.c (! isnan (s2.c))) > 100000000*eps) ...
         || any (abs (s1.e (! isnan (s1.e)) - s2.e (! isnan (s2.e))) > 100000000*eps) ...
         || any (abs (s1.b (! isnan (s1.b)) - s2.b (! isnan (s2.b))) > 100000000*eps) ...
         || ! strcmp(class(s1.c), class(s2.c))
         || ! strcmp(class(s1.e), class(s2.e))
         || ! strcmp(class(s1.b), class(s2.b)))
        x   = ts(j).data;
        fts = ts(j).tests(i);
        fr1 = s1;
        fr2 = s2;
        return;
      endif
    endfor
  endfor
endfunction
