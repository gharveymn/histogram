function [counts,edges,binidx] = histcounts (data, varargin)
  
  MAXBINS = 65536;

  validateattributes (data, {"numeric", "logical"}, {"real"}, "histcounts", 1);
  isedgestrans = false;
  
  if (nargin < 3)
    ins = [];
    if (nargin == 2 && ! isscalar (varargin{1}))
      currin = varargin{1};
      validateattributes (currin, {"numeric", "logical"}, {"vector",...
                         "nonempty", "real", "nondecreasing"}, "histcounts", 2);
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
                            "positive"}, "histcounts", 2);
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
    counts = histcountsmex (data, edges);
  else
    [counts, binidx] = histcountsmex (data, edges);
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

function edges = bmselect (bm, datacol, datamin, datamax, hardlim)
  switch (ins)
    case "auto"
      edges = bmsauto (datacol, datamin, datamax, hardlim);
    case "scott"
      edges = bmscott (datacol, datamin, datamax, hardlim);
    case "fd"
      edges = bmfd (datacol, datamin, datamax, hardlim);
    case "integers"
      edges = bmintegers (datacol, datamin, datamax, hardlim, MAXBINS);
    case "sqrt"
      edges = bmsqrt (datacol, datamin, datamax, hardlim);
    case "sturges"
      edges = bmsturges (datacol, datamin, datamax, hardlim);
  endswitch
endfunction
