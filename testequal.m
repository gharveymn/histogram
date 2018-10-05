function [x, fts, fr1, fr2] = testequal(ts,r1,r2)
  x   = [];
  fts = [];
  fr1 = [];
  fr2 = [];
  for i = 1:numel(r1)
    retval = teststructequal(r1{i},r2{i});
    if (retval != 0)
        x   = ts{i}.data;
        fts = ts{i}.tests(retval);
        fr1 = r1{i}(retval,:);
        fr2 = r2{i}(retval,:);
    endif
  endfor
endfunction
