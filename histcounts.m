function [counts,edges,binidx] = histcounts (data, varargin)
  
  validateattributes (data, {"numeric", "logical"}, {"real"}, "histcounts", 1);
  areedgestrans = false;
  
  if (nargin < 3)
    opts = [];
    if (nargin == 2 && ! isscalar (varargin{1}))
      currin = varargin{1};
      validateattributes (currin, {"numeric", "logical"}, {"vector",...
                         "nonempty", "real", "nondecreasing"}, "histcounts", 2);
      if (iscolumn (currin))
        edges = currin.';
        areedgestrans = true;
      else
        edges = currin;
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
      
      if (nargin == 1)
        edges = bmauto (datacol, datamin, datamax, false);
      else
        currin = varargin{1};
        validateattributes (currin, {"numeric", "logical"}, {"integer",...
                            "positive"}, "histcounts", 2);
        datarange = datamax - datamin;
        numbins   = double(currin);
        edges = pickbins (datamin, datamax, numbins, datarange/numbins);
      endif
    endif
    
  else
    
    
    
  endif
  
endfunction
