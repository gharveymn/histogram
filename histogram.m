classdef histogram < handle
	
  properties
    Annotation         = [];
    BeingDeleted       = "off";
    BinCounts          = [];         ## 0
    BinCountsMode      = [];         ## "auto"
    BinEdges           = [];         ## [0 1]
    BinLimits          = [];         ## [0 1]
    BinLimitsMode      = [];         ## "auto"
    BinMethod          = [];         ## "auto"
    BinWidth           = [];         ## 1
    BusyAction         = "queue";
    ButtonDownFcn      = "";
    Children           = [];
    CreateFcn          = "";
    Data               = [];
    DeleteFcn          = "";
    DisplayName        = "";
    DisplayStyle       = "bar";
    EdgeAlpha          = 1;
    EdgeColor          = [0 0 0];
    FaceAlpha          = 0.6;
    FaceColor          = "auto";
    HandleVisibility   = "on";
    HitTest            = "on";
    Interruptible      = "on";
    LineStyle          = "-";
    LineWidth          = 0.5;
    Normalization      = "";         ## "count"
    NumBins            = [];         ## 1
    Orientation        = "";         ## "vertical"
    Parent             = [];
    PickableParts      = "visible";
    Selected           = "off";
    SelectionHighlight = "on";
    Tag                = "";
    Type               = "histogram";
    UIContextMenu;
    UserData           = [];
    Values             = 0;          ## only here for MATLAB compat
    Visible            = "on";
  endproperties
  
  properties (Hidden, SetAccess = private)
    AutoColor = [];
  endproperties
  
  methods (Access = public)
    
    function h = histogram (varargin)
      
      hax = __plt_get_axis_arg__ ("histogram", varargin{:});
      innum = 1 + (! isempty (hax));
      if (innum <= nargin && ~ischar (varargin{innum}))
        datavarname = inputname (innum);
      else
        datavarname = "";
      endif
      
      dispatch = parseinput (h, varargin, nargin, innum);
      
      if (dispatch)
        ## remove error with the implementation of categorical
        error (["histogram: the 'categorical' class has not yet been "...
                "implemented in octave, cannot continue"]);
        if (isempty (hax))
          h = histogram (categorical([]), varargin{:});
        else
          h = histogram (hax, categorical([]), varargin{2:nargin});
        endif
      endif
      
      hax = newplot (hax);
      
      # determine auto coloring
      clrord = get (hax, 'ColorOrder');
      clridx = [get(hax, 'ColorOrderIndex'), get(hax, 'LineStyleOrderIndex')];
      clrsiz = size (clrord, 1);
      h.AutoColor = clrord(mod (clridx(1) - 1, clrsiz) + 1, :);
      
      x = h.Data;
      if (isempty (x) && ! isempty (h.BinCounts))
        x = h.BinEdges;
      endif
      
      if (strcmp (h.Orientation, "horizontal"))
        ## FIXME: should error here if hax is class 'PolarAxes', whenever 
        ##        that is implemented
        
        ## FIXME: need configureAxes here whenever implemented (not vital 
        ##        though since it doesn't do much here)
      else
        ## FIXME: configureAxes here
      endif
      
      if (isempty (h.BinCounts))
        numnvargs = 0;
        nvargs = cell(8,1);
        if (! isempty (h.Normalization))
          nvargs(1) = "Normalization";
          nvargs(2) = h.Normalization;
          numnvargs += 2;
        endif
        
        if(! isempty (h.BinMethod))
          nvargs(numnvargs+1) = "BinMethod";
          nvargs(numnvargs+2) = h.BinMethod;
          numnvargs += 2;
        endif
        
        if (! isempty (h.BinEdges))
          h.BinCounts = histcounts (h.Data, h.BinEdges, nvargs{1:numnvargs});
        else
          if (! isempty (h.BinLimits))
            nvargs(numnvargs+1) = "BinLimits";
            nvargs(numnvargs+2) = h.BinLimist;
            numnvargs += 2;
          endif
          
          if (! isempty (h.BinWidth))
            nvargs(numnvargs+1) = "BinWidth";
            nvargs(numnvargs+2) = h.BinWidth;
            numnvargs += 2;
          endif
          
          if(! isempty (h.NumBins))
            [h.BinCounts, h.BinEdges] = histcounts (h.Data, h.NumBins,...
                                                    nvargs{1:numnvargs});
          else
            [h.BinCounts, h.BinEdges] = histcounts (h.Data,...
                                                    nvargs{1:numnvargs});
          endif
        endif
      endif
      
      xbars = [h.BinEdges(1:end-1)
               h.BinEdges(1:end-1)
               h.BinEdges(2:end)
               h.BinEdges(2:end)];
      ybars = [zeros(1,numel(h.BinCounts))
               h.BinCounts
               h.BinCounts
               zeros(1,numel(h.BinCounts))];
      
      patch (hax, xbars, ybars, "b");
      
      h.Parent = hax;
      
    endfunction
    
  endmethods
  
  methods (Access = private)
    function dispatch = parseinput (h, args, nargs, innum);
      
      # pi = struct ("Data", [], "NumBins", [], "BinEdges", [],...
      #             "BinLimits", [], "BinWidth", [], "Normalization", "",...
      #             "BinMethod", "auto", "BinCounts", [], "Orientation", "");
      
      hax = __plt_get_axis_arg__ ("histogram", args{:});
      innum = 1 + (! isempty (hax));
      
      foundonlynamevals = true;
      dispatch = false;
      
      ## first input
      if (innum <= nargs)
        currin = args{innum};
        if(! ischar (currin))
          validateattributes (currin, {"numeric", "logical", "datetime",...
                              "duration"}, {"real"}, "histogram", "x", innum);
          h.Data = currin;
          foundonlynamevals = false;
          innum += 1;
        endif
        
        ## second input
        if (innum <= nargs)
          currin = args{innum};
          if (! ischar (currin))
            if (isscalar (currin))
              validateattributes (currin, {"numeric", "logical"}, {"integer",...
                                  "positive"}, "histogram", "nbins", innum);
              h.NumBins   = currin;
              h.BinMethod = "";
            else
              ## unimplemented functions here won't run until datetime and/or
              ## duration are implemented
              if (isa (h.Data, "datetime"))  
                if (! isa (currin, "datetime") || ! isvector (currin)...
                    || isempty (currin))
                    error ("histogram: invalid edges for datetime input");
                elseif (! issorted (currin) || any (isnat (currin)))
                    error ("histogram: unsorted eges for datatime input");
                endif    
              elseif (isa(h.Data, "duration"))
                if (! isa (currin, "duration") || ! isvector (currin)...
                    || isempty (currin))
                    error ("histogram: invalid edges for duration input");
                elseif (! issorted (currin) || any (isnan (currin)))
                    error ("histogram: unsorted eges for duration input");
                endif
              else
                validateattributes (currin, {"numeric", "logical"}, {"vector",...
                                    "nonempty", "real", "nondecreasing"}, ...
                                    "histogram", "edges", innum);
              endif
              h.BinEdges  = currin;
              h.BinMethod = "";
            endif
            innum += 1;
          endif
          
          ## parse name-val pairs
          if (rem (nargs-innum+1, 2) != 0)
            error ("histogram: each name argument must have a value pair");
          endif
          
          ## does not include Annocation, BeingDeleted, Children, Type, Values
          validnames = {"BinCounts",
                        "BinCountsMode",
                        "BinEdges",
                        "BinLimits",
                        "BinLimitsMode",
                        "BinMethod",
                        "BinWidth",
                        "BusyAction",
                        "ButtonDownFcn",
                        "CreateFcn",
                        "Data",
                        "DeleteFcn",
                        "DisplayName",
                        "DisplayStyle",
                        "EdgeAlpha",
                        "EdgeColor",
                        "FaceAlpha",
                        "FaceColor",
                        "HandleVisibility",
                        "HitTest",
                        "Interruptible",
                        "LineStyle",
                        "LineWidth",
                        "Normalization",
                        "NumBins",
                        "Orientation",
                        "Parent",
                        "PickableParts",
                        "Selected",
                        "SelectionHighlight",
                        "Tag",
                        "UIContextMenu",
                        "UserData",
                        "Visible"};
          
          for j = innum:2:nargs
            name  = validatestring (args{j}, validnames);
            value = args{j+1};
            
            switch (name)
              case "Data"
                validateattributes (value, {"numeric", "logical", "datetime",...
                              "duration"}, {"real"}, "histogram", "Data", j+1);
                h.Data = value;
              case "NumBins"
                validateattributes (value, {"numeric", "logical"}, {"scalar",...
                                    "integer", "positive"}, "histogram",...
                                    "NumBins", j+1);
                if (! isempty (h.BinEdges))
                  error ("histogram: invalid mixture of bin arguments");
                endif
                h.NumBins   = value;
                h.BinMethod = "";
                h.BinWidth  = [];
              case "BinEdges"
                ## unimplemented functions here won't run until datetime and/or
                ## duration are implemented
                if (isa (value, "datetime"))  
                  if (! isvector (value))
                      error ("histogram: invalid edges for datetime input");
                  elseif (! issorted (value) || any (isnat (value)))
                      error ("histogram: unsorted eges for datatime input");
                  endif    
                elseif (isa(value, "duration"))
                  if (! isvector (value))
                      error ("histogram: invalid edges for duration input");
                  elseif (! issorted (value) || any (isnan (value)))
                      error ("histogram: unsorted eges for duration input");
                  endif
                else
                  validateattributes (value, {"numeric", "logical"}, {"vector",...
                                      "nonempty", "real", "nondecreasing"}, ...
                                      "histogram", "BinEdges", j+1);
                endif
                
                if (length (value) < 2)
                  error("histogram: must have at least 2 bin edges.");
                endif
                
                h.BinEdges  = value;
                h.BinMethod = "";
                h.BinWidth  = [];
                h.NumBins   = [];
                h.BinLimits = [];
              case "BinWidth"
                if (isa (value, "duration"))  
                  if (! isscalar (value) && isfinite(value))
                        error (["histogram: invalid bin width for duration "...
                               "input"]);
                  endif    
                elseif (isa(value, "calendarDuration"))
                  if (! isscalar (value) && isfinite(value))
                        error (["histogram: invalid bin width for "...
                                "calendarDuration input"]);
                  endif
                  [caly,calm,cald,calt] = split (value,...
                                                 {'year','month','day','time'});
                  if (caly < 0 || calm < 0 || cald < 0 || calt < 0) || ...
                          (caly == 0 && calm == 0 && cald == 0 && calt == 0)
                      error (["histogram: invalid bin width for "...
                                "calendarDuration input"]);
                  end
                else
                  validateattributes (value, {"numeric", "logical"}, {"scalar",...
                                      "real", "positive", "finite"}, ...
                                      "histogram", "BinWidth", j+1);
                endif
                if (! isempty (h.BinEdges))
                  error ("histogram: invalid mixture of bin arguments");
                endif
                h.BinWidth  = value;
                h.BinMethod = "";
                h.NumBins   = [];
              case "BinLimits"
                if (isa (value, "datetime") || isa (value, "duration"))
                  if (! (numel (value) == 2 && issorted (value)...
                         && all (isfinite (value))))
                    error("histogram: invalid datetime or duration bin limits");     
                  endif
                else
                  validateattributes (value, {"numeric", "logical"}, {"numel",...
                                      2, "vector", "real", "nondecreasing",...
                                      "finite"}, "histogram", "BinLimits", j+1);                    
                endif
                if (! isempty (h.BinEdges))
                  error ("histogram: invalid mixture of bin arguments");
                endif
                h.BinLimits = value;
              case "Normalization"
                h.Normalization = validatestring (value, {"count",...
                                                  "countdensity", "cumcount",...
                                                  "probability", "pdf",...
                                                  "cdf"}, "histogram",...
                                                  "Normalization", j+1);
              case "BinMethod"
                h.BinMethod = validatestring (value, {"auto", "scott", "fd",...
                                              "integers", "sturges", "sqrt",...
                                              "century", "decade", "year",...
                                              "quarter", "moneth", "week",...
                                              "day", "hour", "minute",...
                                              "second"}, "histogram",...
                                              "BinMethod", j+1);
                if (! isempty (h.BinEdges))
                  error ("histogram: invalid mixture of bin arguments");
                endif
                h.BinWidth = [];
                h.NumBins  = [];
              case "BinCounts"
                validateattributes (value, {"numeric", "logical"}, {"real",...
                                    "vector", "nonnegative", "finite"},...
                                    "histogram", "BinCounts", j+1);
                h.BinCounts = reshape(value, 1, []);                    
              case "BinLimitsMode"
                h.BinLimitsMode = validatestring(value, {"auto", "manual"});
              case "BinCountsMode"
                h.BinCountsMode = validatestring(value, {"auto", "manual"});
              case "Categories"
                if (foundonlynamevals)
                  dispatch = true;
                  return;
                else
                  error ("histogram: categories unsupported with this syntax");
                endif
              case "Orientation"
                h.Orientation = validatestring (value, {"vertical",...
                                                "horizontal"}, "histogram",...
                                                "Orientation", j+1);
              otherwise
                h.(name) = value;
            endswitch
            
          endfor
          
          ## consistency checks
          
          if (! isempty (h.BinLimitsMode))
            if (strcmp (h.BinLimitsMode, "auto"))
              if (! isempty (h.BinEdges) || ! isempty (h.BinLimits))
                error("histogram: non-empty bin limits or edges in auto mode");
              endif
            else
              if (isempty (h.BinEdges) && isempty (h.BinLimits))
                error("histogram: empty bin limits or edges in manual mode");
              endif
            endif
          endif
          
          if (! isempty (h.BinCountsMode))
            if (strcmp (h.BinCountsMode, "auto"))
              if (! isempty (h.BinCounts))
                error("histogram: non-empty bin counts in auto mode");
              endif
            else
              if (isempty (h.BinCounts))
                error("histogram: empty bin counts in auto mode");
              endif
            endif
          endif
          
          if (! isempty (h.BinCounts))
            if (length (h.BinEdges) != length (h.BinCounts)+1)
              error(["histogram: the length of BinEdges must equal one "...
                     "more than BinCounts"]);
            endif
            
            if (! isempty (h.Data))
              error("histogram: cannot mix data and specified bin counts");
            endif
            
            if ((isa (h.BinEdges, "datetime") || isa (h.BinEdges, "duration"))...
                && ismember (h.Normalization, {"countdensity", "pdf"}))
                error ("histogram: datetime or duration normalization error");
            endif
          else
            if (isa (h.Data, "datetime") || isa (h.Data, "duration"))
              if (! isempty (h.BinEdges) && ! isequal (class (h.Data),...
                  class (h.BinEdges)))
                error (["histogram: the datetime bin edges class '"...
                        class(h.BinEdges) "' does not match the Data class '"...
                        class(h.Data) "'"]);
              elseif (! isempty (h.BinLimits) && ! isequal (class (h.Data),...
                      class (h.BinLimits)))
                error (["histogram: the datetime bin edges class '"...
                        class(h.BinLimits) "' does not match the Data class '"...
                        class(h.Data) "'"]);  
              elseif (ismember (h.Normalization, {"countdensity", "pdf"}))
                error (["histogram: invalid datetime or duration "
                        "normalization 'countdensity' or 'pdf'"]);
              elseif (ismember (h.BinMethod, {"integers"}))
                error (["histogram: invalid datetime or duration BinMethod "...
                        "'integers'"]);
              endif
              
              if (isa (h.Data, "datetime"))
                if (! isempty (h.BinWidth) && (isnumeric (h.BinWidth)...
                    || islogical (h.BinWidth)))
                  error ("histogram: invalid BinWidth for datetime Data");
                endif
              else
                if (ismember (h.BinMethod, {"century", "decade", "quarter",...
                    "month", "week"}))
                  error ("histogram: invalid BinMethod for duration Data");
                endif
                if (! isempty (h.BinWidth) && ! isa (h.BinWidth, "duration"))
                  error (["histogram: BinWidth must be of class of "...
                          "'duration' for duration Data"]);
                endif
              endif
              
            else
              
              if (isa (h.BinEdges, "datetime") || isa (h.BinEdges, "duration"))
                error (["histogram: unexpected '" class(h.BinEdges)... 
                        "' BinEdges for numeric Data"]);
              elseif (isa (h.BinWidth, "duration") ...
                      || isa (h.BinWidth, "calendarDuration"))
                error (["histogram: unexpected '" class(h.BinWidth)... 
                        "' BinWidth for numeric Data"]);  
              elseif (isa (h.BinLimits, "datetime")...
                      || isa (h.BinLimits, "duration"))
                error (["histogram: unexpected '" class(h.BinLimits)... 
                        "' BinLimits for numeric Data"]);
              elseif (ismember (h.BinMethod, {"century", "decade", "year",...
                      "quarter", "month", "week", "day", "hour", "minute",...
                      "second"}))
                error (["histogram: invalid numeric BinMethod"]);
              endif
            endif
          endif
        endif
      endif
    endfunction
  endmethods
endclassdef
