## Copyright (C) 2018-2018 Gene Harvey
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {[@var{counts}, @var{edges}, @var{binidx}] =} histcounts (@var{A})
## @deftypefnx {} {[@var{counts}, @var{edges}, @var{binidx}] =} histcounts (@var{A}, @var{nbins})
## @deftypefnx {} {[@var{counts}, @var{edges}, @var{binidx}] =} histcounts (@var{A}, @var{edges})
## @deftypefnx {} {[@var{counts}, @var{edges}, @var{binidx}] =} histcounts (@var{A}, @dots{}, @var{prop}, @var{val}, @dots{})
## Calculate histogram binning counts.
##
## In the simplest case, one may get the binning counts with 
##
## @example
## @var{counts} = histcounts (@var{A})
## @end example
##
## @noindent
## where @var{A} can by any numeric N-dimensional array, and @var{counts} is a 
## row vector.
##
## The edges used to calculate bin counts are returned with 
##
## @example
## [@var{counts}, @var{edges}] = histcounts (@var{A});
## @end example
##
## @noindent
## where @var{edges} is a row vector of length @code{numel (@var{counts}) + 1}.
##
## The bin index of each element in @var{A} is returned with 
##
## @example
## [@var{counts}, @var{edges}, @var{binidx}] = histcounts (@var{A});
## @end example
##
## @noindent
## where @var{binidx} has the same dimensions as @var{A}.
##
## The number of bins may be specified with a positive numeric or logical 
## integer scalar.  The syntax for this is 
##
## @example
## @var{counts} = histcounts (@var{A}, @var{nbins});
## @end example
##
## The bin edges may be specified by a non-decreasing numeric or logical 
## vector.  Accomplish this with 
##
## @example
## @var{counts} = histcounts (@var{A}, @var{edges});
## @end example
## 
## @noindent
## Specifying the bin edges is incompatible with all properties except 
## the @qcode{"Normalization"} property.
##
## Properties:
##
## @table @asis
## @item "NumBins"
## The number of bins to use.  This must be a positive numeric or logical 
## integer scalar.
##
## @item "BinEdges"
## The edges for which to bin the input @var{A}.  This must be a non-decreasing 
## numeric or logical vector. This property is incompatible with all properties 
## except the @qcode{"Normalization"} property.
##
## @item "BinWidth"
## The width of bins to use.  A limit of 65536 bins is imposed to prevent 
## creation of too many bins.  If this limit is hit, then the bin width
## is increased to accommodate for this.
##
## @item "BinLimits"
## The minimum and maximum of created bins.  Counts will be taken only from 
## elements in @var{A} between the limits (inclusive).  This must be a two-
## element non-decreasing numeric or logical vector.
##
## @item "Normalization"
## The normalization algorithm applied to @var{counts} before return.  The 
## available algorithms are:
##
## @table @asis
## @item "count"
## [default] The number of counts in each bin using standard histogram 
## technique.
##
## @item "probability"
## The probability of an observation being in each bin. Each bin is divided by 
## the number of elements in @var{A}.  Note that the number of elements 
## includes elements that may not have been counted because they were outside 
## of the range of bin edges.
##
## @item "countdensity"
## The density of counts in each bin.  The counts in each bin are divided by 
## the bin width (element-wise).
##
## @item "pdf"
## An estimate of the probability density function.  This is the number counts 
## in each bin divided by the total number of observations divided by the 
## bin width (element-wise).
##
## @item "cumcount"
## A cumulative sum given by the addition of each bin and all previous bins.
##
## @item "cdf"
## An estimate of the cumulative density function.  This is a cumulative sum 
## of the bins divided by the number of elements in @var{A}.  Note that the 
## number of elements includes elements that may not have been counted because 
## they were outside of the range of bin edges.
##
## @end table
##
## @item "BinMethod"
## The method used to automatically calculate edges for the binning of @var{A}.
## The available methods are:
##
## @table @asis
## @item "auto"
## [default] Use an automatic algorithm for calculating the bin edges.  If 
## @var{A} is integer data, then this will most likely use the 
## @qcode{"integers"} method.  Otherwise it will use the @qcode{"scott"} 
## method.
##
## @item "scott"
## Use Scott's rule to calculate the bin edges.  The edges are calculated with
## @code{7/2 * std (@var{A}(:)) * numel (@var{A}) ^ (-1 / 3)}.  This is ideal
## for data which is close to Gaussian.
##
## @item "fd"
## Use the Freedman-Diaconis rule. The edges are calculated with 
## @code{2 * iqr (@var{A}(:)) * numel (@var{A}) ^ (-1 / 3)}.  This method is 
## ideal for data with outliers.
##
## @item "integers"
## Create an individual bin for each point of integer data.  The bins are 
## offset by 1/2 whole number.  This is the case unless the range of @var{A} 
## is greater than 65536, in which case lerger bin widths are used.
##
## @item "sturges"
## Use Sturges' rule to calculate bin edges.  The edges are calcuated with 
## @code{ceil (1 + log2 (numel (@var{A})))}.
##
## @item "sqrt"
## Use the square root rule to calculate bin edges.  The edges are calculated 
## with @code{ceil (sqrt (numel (@var{A})))}.
## @end table
## @end table
##
## @seealso{histc, hist}
## @end deftypefn

function [counts,edges,binidx] = histcounts (data, varargin)
 
  MAXBINS = 65536;

  validateattributes (data, {"numeric", "logical"}, {"real"}, "histcounts");
  isedgestrans = false;
  
  s  = [];
  si = [];
  
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
          [edges,s,si] = bmselect (ins.BinMethod, datacol, datamin, datamax, false);
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
        edges = [datamin + (0:numbins-1) .* ((datamax - datamin) / numbins),...
                 datamax];
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
  
  if (nargout <= 2)
    counts = histcountsoct (data, edges, ! isempty(s));
  else
    if (! isempty (s))
      s = reshape(s, size(data));
      [counts, binidx] = histcountsoct (s, edges, true);
      binidx(si) = binidx;
    else
      [counts, binidx] = histcountsoct (data, edges, false);
    endif
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
                            "integer", "positive"}, "histcounts", "NumBins", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.NumBins   =  value;
        opts.BinMethod = [];
        opts.BinWidth  = [];
      case "BinEdges"
        validateattributes (value, {"numeric", "logical"}, {"vector", "real",...
                            "nondecreasing"}, "histcounts", "BinEdges", j + 2);
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
                            "positive", "finite"}, "histcounts", "BinWidth", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinWidth  = value;
        opts.BinMethod = [];
        opts.NumBins   = [];
      case "BinLimits"
        validateattributes (value, {"numeric", "logical"}, {"numel", 2,...
                            "vector", "real", "nondecreasing", "finite"},...
                            "histcounts", "BinLimits", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinLimits = value;
      case "Normalization"
        opts.Normalization = validatestring (value, {"count", "countdensity",...
                                             "cumcount", "probability",...
                                             "pdf", "cdf"}, "histcounts",...
                                             "Normalization", j + 2);
      otherwise ## aka "BinMethod"
        opts.BinMethod = validatestring (value, {"auto", "scott", "fd",...
                                         "integers", "sturges", "sqrt"},...
                                         "histcounts", "BinMethod", j + 2);
        if (! isempty (opts.BinEdges))
          error ("histcounts: invalid mixture of binning inputs");
        endif
        opts.BinWidth = [];
        opts.NumBins  = [];
    endswitch
  endfor
endfunction

function [edges, s, si] = bmselect (bm, d, dl, dh, haslim)
  s  = [];
  si = [];
  switch (bm)
    case "auto"
      edges = bmauto (d, dl, dh, haslim);
    case "scott"
      edges = bmscott (d, dl, dh, haslim);
    case "fd"
      [edges,s,si] = bmfd (d, dl, dh, haslim);
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
    edges = bmintegers (d, dl, dh, haslim, MAXBINS);
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

function [edges, s, si] = bmfd (d, dl, dh, haslim)
  n = numel (d);
  if (n > 1)
    [s, si] = sort(d);
    dr = s(end) - s(1);
    iq = max (sortediqr (s(:)), double (dr) / 10);
    bw = 2 * iq * n ^ (-1 / 3);
    
    if (! haslim)
      edges = pickbins (dl, dh, [], bw);
    else
      edges = pickbinsbl (s(1), s(end), dl, dh, bw);
    endif
    
  else
    dr = 0;
    bw = 1;
    
    s = [];
    si = [];
    
    if (! haslim)
      edges = pickbins (dl, dh, [], bw);
    else
      edges = pickbinsbl (d, d, dl, dh, bw);
    endif
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
        
        numbins = innumbins;
        right = min (max ( left + numbins .* bw, dh), realmax (class (dh))); 
        
      endif
      
    else
      if (isempty (innumbins))
        innumbins = 1;
      endif
      
      br    = max (1, ceil (innumbins * eps (ds)));
      left  = floor (2 * (dl - br ./ 4)) / 2;
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

function iq = sortediqr(x)
  n = sum(! isnan (x(:)));
  iq = sortedprctile(x, 75, n) - sortedprctile(x, 25, n);
endfunction

function v = sortedprctile(x, p, n)
  r = (p/100)*n;
  k = floor(r+0.5);
  kp1 = k+1;
  r = r-k;
  
  k   = max (k, 1);
  kp1 = min (kp1, n);
  
  if(k != kp1)
    v = (0.5 + r) .* x(kp1) + (0.5 - r) .* x(k);
  else
    v = x(k);
  end
endfunction

function n = MAXBINS
  n = 65536;
endfunction
