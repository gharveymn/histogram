function [counts,edges,binidx] = histcounts (data, varargin)
 
  MAXBINS = 65536;

  validateattributes (data, {"numeric", "logical"}, {"real"}, "histcounts");
  isedgestrans = false;
  
  if (nargin < 3)
    ins = [];
    if (nargin == 2 && ! isscalar (varargin{1}))
      currin = varargin{1};
      validateattributes (currin, {"numeric", "logical"}, {"vector",...
                         "nonempty", "real", "nondecreasing"}, "histcounts");
      if (iscolumn (currin))
        edges = currin.';
        isedgestrans = true;
      else
        edges = currin;
      endif
    else
      if (! isfloat (data))
        datacol = data(:);
        datamin = double (min (datacol));
        datamax = double (max (datacol));        
      else
        datacol = data(isfinite(data));
        datamin = min (datacol);
        datamax = max (datacol);
      endif
      
      if (nargin == 1)
        edges = bmauto (datacol, datamin, datamax, false);
      else
        currin = varargin{1};
        validateattributes (currin, {"numeric", "logical"}, {"integer",...
                            "positive"}, "histcounts");
        datarange = datamax - datamin;
        numbins   = double (currin);
        edges = pickbins (datamin, datamax, numbins, datarange/numbins);
      endif
    endif
  else
    ins = parseinput(varargin);

    if (isempty (ins.BinLimits))
      if (! isempty (ins.BinEdges))
        if (iscolumn (ins.BinEdges))
          edges = ins.BinEdges.';
          isedgestrans = true;
        else
          edges = ins.BinEdges;
        endif
      else
        if (! isfloat (data))
          datacol = data(:);
          datamin = double(min(datacol));
          datamax = double(max(datacol));        
        else
          datacol = data(isfinite(data));
          datamin = min (datacol);
          datamax = max (datacol);       
        endif

        if (! isempty (ins.NumBins))
          datarange = datamax - datamin;
          numbins   = double (ins.NumBins);
          edges = pickbins (datamin, datamax, numbins, datarange/numbins);
        elseif (! isempty (ins.BinWidth))
          if (! isfloat (ins.BinWidth))
            ins.BinWidth = double (ins.BinWidth);
          endif
          datarange = datamax - datamin;
          if (! isempty (datamin))
            left = ins.BinWidth * floor (datamin / ins.BinWidth);
            numbins = max (1, ceil ((datamax - left) ./ ins.BinWidth));
            if (numbins > MAXBINS)
              ins.BinWidth = datarange/(MAXBINS-1);
              left = ins.BinWidth * floor (datamin / ins.BinWidth);
              if (datamax <= left + (numbins - 1) * ins.BinWidth)
                ins.BinWidth = datarange / MAXBINS;
                left = datamin;
              endif
            endif
            edges = left + (0:numbins) .* ins.BinWidth;
          else
            edges = cast ([0, ins.BinWidth], class (datarange));
          endif
        else
          edges = bmselect (ins.BinMethod, datacol, datamin, datamax, false);
        endif
      endif
    else
      if (! isfloat (ins.BinLimits))
        datamin = double (ins.BinLimits(1));
        datamax = double (ins.BinLimits(2));
      else
        datamin = ins.BinLimits(1);
        datamax = ins.BinLimits(2);
      endif
      if (! isempty (ins.NumBins))
        numbins = double (ins.NumBins);
        edges = [datamin + (0:numbins-1) .* ((datamax - datamin) / numbins), datamax];
      elseif (! isempty (ins.BinWidth))
        if (! isfloat (ins.BinWidth))
          ins.BinWidth = double (ins.BinWidth);
        endif
        ins.BinWidth = max (ins.BinWidth, (datamax - datamin) / MAXBINS);
        edges = datamin:ins.BinWidth:datamax;
        if (edges(end) < datamax || isscalar (edges))
          edges = [edges, datamax];
        endif
      else
        datacol = data(data >= datamin & data <= datamax);
        edges = bmselect (ins.BinMethod, datacol, datamin, datamax, true);
      endif
    endif
  endif

  edges = full(edges);

  ## TODO histcountsmex in octfile
  if (nargout <= 2)
    counts = histcountsnaive (data, edges);
  else
    [counts, binidx] = histcountsnaive (data, edges);
  endif

  if (! isempty (ins))
    switch ins.Normalization
      case "countdensity"
        counts = counts ./ double (diff (edges));
      case "cumcount"
        counts = cumsum (counts);
      case "probability"
        counts = counts / numel (data);
      case "pdf"
      counts = counts / numel (data) ./ double (diff (edges));
      case "cdf"
        counts = cumsum (counts / numel (data));
    endswitch
  endif

  if (nargin > 1 && isedgestrans)
    edges = edges.';
  endif
endfunction

function opts = parseinput (input)
  
  opts = struct ("NumBins", [], "BinEdges", [], "BinLimits", [],...
                 "BinWidth", [], "Normalization", "count",...
                 "BinMethod", "auto");
  
  idx = 1;
  currin = input{idx};
  if (isnumeric (currin) || islogical (currin))
    if (isscalar (currin))
      validateattributes (currin, {"numeric", "logical"}, {"integer",...
                          "positive"}, "histcounts", idx + 1);
      opts.NumBins   = currin;
      opts.BinMethod = [];
    else
      validateattributes (currin, {"numeric", "logical"}, {"vector", ...
                          "nonempty", "real", "nondecreasing"},...
                          "histcounts", idx + 1);
      opts.NumEdges  = currin;
      opts.BinMethod = [];
    endif
    idx += 1;
  endif
  
  ## parse name-val pairs
  if (rem (length(input)-idx+1, 2) != 0)
    error ("histcount: each name argument must have a value pair");
  endif
  
  validnames = {"NumBins", "BinEdges", "BinWidth", "BinLimits",...
                "Normalization", "BinMethod"};
  
  for j = idx:2:length(input)
    name  = validatestring (input{j}, validnames, "histcounts", j + 1);
    value = input{j+1};
    switch (name)
      case "NumBins"
        validateattributes (value, {"numeric", "logical"}, {"scalar",...
                            "integer", "positive"}, "NumBins", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.NumBins   =  value;
        opts.BinMethod = [];
        opts.BinWidth  = [];
      case "BinEdges"
        validateattributes (value, {"numeric", "logical"}, {"vector", "real", 
                            "nondecreasing"}, "histcounts", j + 2);
        if (length (value) < 2)
          error("histcounts: must have at least 2 bin edges.");
        endif
        opts.BinEdges  = value;
        opts.BinMethod = [];
        opts.NumBins   = [];
        opts.BinWidth  = [];
        opts.BinLimits = [];
      case "BinWidth"
        validateattributes (value, {"numeric", "logical"}, {"scalar", "real",...
                            "positive", "finite"}, "BinWidth", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinWidth  = value;
        opts.BinMethod = [];
        opts.NumBins   = [];
      case "BinLimits"
        validateattributes (value, {"numeric", "logical"}, {"numel", 2,...
                            "vector", "real", "nondecreasing", "finite"},...
                            "histcounts", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinLimits = value;
      case "Normalization"
        opts.Normalization = validatestring (value, {"count", "countdensity",...
                                             "cumcount", "probability",...
                                             "pdf", "cdf"}, "histcounts",...
                                             j + 2);
      otherwise ## aka "BinMethod"
        opts.BinMethod = validatestring (value, {"auto", "scott", "fd",...
                                         "integers", "sturges", "sqrt"},...
                                         "histcounts", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinWidth = [];
        opts.NumBins  = [];
    endswitch
  endfor
endfunction

function edges = bmselect (bm, d, dl, dh, haslim)
  switch (bm)
    case "auto"
      edges = bmauto (d, dl, dh, haslim);
    case "scott"
      edges = bmscott (d, dl, dh, haslim);
    case "fd"
      edges = bmfd (d, dl, dh, haslim);
    case "integers"
      edges = bmintegers (d, dl, dh, haslim, MAXBINS);
    case "sqrt"
      edges = bmsqrt (d, dl, dh, haslim);
    case "sturges"
      edges = bmsturges (d, dl, dh, haslim);
  endswitch
endfunction

function edges = bmauto (d, dl, dh, haslim)
  dr = dh - dl;
  if (! isempty (d) && (isinteger (d) || islogical (d)...
      || isequal (round (d), d)) && dr <= 50 ...
      && dh <= flintmax (class (dh)) / 2 && dl >= -flintmax( class (dl)) / 2)
    edges = bmintegers (d, dl, dh, haslim, 65536); ## refer to MAXBINS
  else
    edges = bmscott (d, dl, dh, haslim);
  endif  
endfunction

function edges = bmscott (d, dl, dh, haslim)
  if (! isfloat (d))
    d = double (d);
  endif
  bw = 3.5 * std (d) / (numel (d) ^ (1/3));
  if (! haslim)
    edges = pickbins (dl, dh, [], bw);
  else
    edges = pickbinsbl (min (d(:)), max (d(:)), dl, dh, bw);
  endif
endfunction

function edges = bmfd (d, dl, dh, haslim)
  n = numel (d);
  dr = max (d(:)) - min (d(:));
  if (n > 1)
    iq = max (iqr (d(:)), double (dr) / 10);
    bw = 2 * iq * n ^ (-1 / 3);
  else
    binwidth = 1;
  endif
  if (! haslim)
    edges = pickbins (dl, dh, [], bw);
  else
    edges = pickbinsbl (min (d(:)), max (d(:)), dl, dh, bw);
  endif
endfunction

function edges = bmintegers (d, dl, dh, haslim, maxbins)
  
  if (! isempty (dh) && (dh > flintmax (class (dh)) / 2 ...
      || dl < -flintmax (class (dl)) / 2))
    error ("histcounts: input out of int range");
  endif
  
  dr = dh - dl;
  if (! isempty (d))
    ds = max (abs (d(:)));
    dr = max (d(:)) - min (d(:));
    if  (dr > maxbins)
      bw = 10 ^ ceil (log10 (dr / maxbins));
    elseif (isfloat (d) && eps (ds) > 1)
      bw = 10 ^ ceil (log10 (eps (ds)));
    else
      bw = 1;
    endif
    if (! haslim)
      dl = bw * round (dl / bw);
      dh = bw * round (dh / bw);
      edges = (floor (dl) - 0.5 * bw):bw:(ceil (dh) + 0.5 * bw);
    else
      dlh = bw * ceil (dl / bw) + 0.5;
      dhl = bw * floor (dh / bw) - 0.5;
      edges = [dl, (dlh:bw:dhl), dh];
    endif
  else
    if (! haslim)
      edges = cast ([-0.5, 0.5], class(dr));
    else
      dlh = ceil (dl) + 0.5;
      dhl = floor (dh) - 0.5;
      edges = [dl, (dlh:bw:dhl), dh];
    endif
  endif
endfunction

function edges = bmsqrt (d, dl, dh, haslim)
  numbins = max (ceil (sqrt (numel (d))), 1);
  if (! haslim)
    bw = (dh-dl)/numbins;
    if (isfinite (bw))
      edges = pickbins (dl, dh, [], bw);
    else
      edges = pickbins (dl, dh, numbins, bw);
    endif
  else
    edges = linspace (dl, dh, numbins + 1);
  endif
endfunction

function edges = bmsturges (d, dl, dh, haslim)
  numbins = max ( ceil (log2 (numel (d)) + 1), 1);
  if (! haslim)
    bw = (dh-dl)/numbins;
    if (isfinite (bw))
      edges = pickbins (dl, dh, [], bw);
    else
      edges = pickbins (dl, dh, numbins, bw);
    endif
  else
    edges = linspace (dl, dh, numbins + 1);
  endif
endfunction

function edges = pickbins (dl, dh, innumbins, inbw)
  if (! isempty (dl))
    ds = max (abs ([dl, dh]));
    dr = dh - dl;
    inbw = max (inbw, eps (ds));
    
    if (dr > max (sqrt (eps (ds)), realmin (class (ds))))
      ord   = 10 .^ floor (log10 (inbw));
      relsz = inbw / ord;
      
      if (isempty (innumbins))
        if (relsz < 1.5)
          bw = 1*ord;
        elseif (relsz < 2.5)
          bw = 2*ord;
        elseif (relsz < 4)
          bw = 3*ord;
        elseif (relsz < 7.5)
          bw = 5*ord;
        else
          bw = 10*ord;
        endif
        
        left = max (min (bw*floor (dl ./ bw), dl), -realmax (class (dh)));
        numbins = max (1, ceil ((dh - left) ./ bw));
        right = min ( max (left + numbins .* bw, dh), realmax (class (dh)));
        
      else
        bw = ord * floor (relsz);
        left = max (min (bw * floor (dl ./ bw), dl), -realmax (class (dl)));
        if (innumbins > 1)
          leftl = (dh - left) / innumbins;
          lefth = (dh - left) / (innumbins - 1);
          leftord = 10 ^ floor (log10 (lefth - leftl));
          bw = leftord * ceil (leftl ./ leftord);
        endif
        
        numbins = innumbins
        right = min (max ( left + numbins .* bw, dh), realmax (class (dh))); 
        
      endif
      
    else
      if (isempty (innumbins))
        innumbins = 1;
      endif
      
      br = max (1, ceil (innumbins * eps (ds)));
      left = floor (2 * (dl - br ./ 4)) / 2;
      right = ceil (2 * (dh + br ./ 4)) / 2;
      
      bw = (right - left) ./ innumbins;
      numbins = innumbins;
      
    endif
    
    if (! isfinite (bw))
      edges = linspace (left, right, numbins + 1);
    else
      edges = [left, left + (1:numbins-1) .* bw, right];
    endif
    
  else
    
    if (! isempty (innumbins))
      edges = cast (0:numbins, class (dl));
    else
      edges = cast ([0,1], class (dl));
    endif
    
  endif
endfunction

function edges = pickbinsbl (dl, dh, llim, hlim, bw)
  ds = max (abs ([dl,dh]));
  dr = dh - dl;
  bw = max (bw, eps (ds));
  
  if (! isempty (dl) && dr > max (sqrt (eps (ds)), realmin (class (ds))))
    numbins = max (ceil ((hlim - llim) / bw), 1);
    edges = linspace (llim, hlim, numbins + 1);
  else
    edges = [llin,hlim];
  endif
  
endfunction

function [counts, binidx] = histcountsnaive (data, edges)
  counts_sz = numel(edges) - 1;
  counts = zeros(1, counts_sz);
  if (nargout == 2)
    binidx = zeros(size(data));
    for i = 1:numel(data)
      if (data(i) <= edges(counts_sz+1))
        for j = counts_sz:-1:1
          if (edges(j) <= data(i))
            counts(j)++;
            binidx(i) = j;
            break;
          endif
        endfor
      endif
    endfor
  else
    for i = 1:numel(data)
      if(data(i) <= edges(counts_sz+1))
        for j = counts_sz:-1:1
          if(edges(j) <= data(i))
            counts(j)++;
            break;
          endif
        endfor
      endif
    endfor
  endif
endfunction
