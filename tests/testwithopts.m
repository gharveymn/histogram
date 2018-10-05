function res = testwithopts(data, testopts)
  res(1,4) = struct('c', [], 'e', [], 'b', []);
  
  numargs = 0;
  args = cell(1,14);
  
  if(~isempty(testopts.NumBinsRaw))
    args{numargs+1} = testopts.NumBinsRaw;
    numargs = numargs+1;
  end
  
  if(~isempty(testopts.BinEdgesRaw))
    args{numargs+1} = testopts.BinEdgesRaw;
    numargs = numargs+1;
  end
  
  if(~isempty(testopts.NumBins))
    args{numargs+1} = 'NumBins';
    args{numargs+2} = testopts.NumBins;
    numargs = numargs+2;
  end
  
  if(~isempty(testopts.BinEdges))
    args{numargs+1} = 'BinEdges';
    args{numargs+2} = testopts.BinEdges;
    numargs = numargs+2;
  end

  if(~isempty(testopts.BinLimits))
    args{numargs+1} = 'BinLimits';
    args{numargs+2} = testopts.BinLimits;
    numargs = numargs+2;
  end
  
  if(~isempty(testopts.BinWidth))
    args{numargs+1} = 'BinWidth';
    args{numargs+2} = testopts.BinWidth;
    numargs = numargs+2;
  end
  
  if(~isempty(testopts.Normalization))
    args{numargs+1} = 'Normalization';
    args{numargs+2} = testopts.Normalization;
    numargs = numargs+2;
  end

  if(~isempty(testopts.BinMethod))
    args{numargs+1} = 'BinMethod';
    args{numargs+2} = testopts.BinMethod;
    numargs = numargs+2;
  end
  
  % do tests with 0, 1, 2, and 3 returns
  histcounts(data, args{1:numargs});
  res(1).c = ans;
  
  [res(2).c] = histcounts(data, args{1:numargs});
  
  [res(3).c, res(3).e] = histcounts(data, args{1:numargs});
  
  [res(4).c, res(4).e, res(4).b] = histcounts(data, args{1:numargs});
end
