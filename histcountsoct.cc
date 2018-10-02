#include <octave/oct.h>
#include <limits>

#define HCO_DEBUG 0

static int hco_CheckEvenSpacing(Matrix in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  double spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz+1; i++)
  {
    if(fabs(spacing - (in_edges(i) - in_edges(i-1))) > std::numeric_limits<double>::epsilon())
    {
      return 0;
    }
  }
  return 1;
}

static int hco_CheckEvenSpacing(FloatMatrix in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  float spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz+1; i++)
  {
    if(fabsf(spacing - (in_edges(i) - in_edges(i-1))) > std::numeric_limits<float>::epsilon())
    {
      return 0;
    }
  }
  return 1;
}

template <class C>
static int hco_CheckEvenSpacing(C in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  typename C::element_type spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz+1; i++)
  {
    if(spacing != (in_edges(i) - in_edges(i-1)))
    {
      return 0;
    }
  }
  return 1;
}

template <class CD, class CE>
static octave_value_list hco_CountNormal(CD in_data, CE in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  Matrix binidx;
  octave_idx_type i, j;
  
  #if HCO_DEBUG
    printf("ran hco_CountNormal\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    binidx = Matrix(in_data.dims());
    binidx.fill(0);
  }
  
  for(i = 0; i < data_sz; i++)
  {
    if((typename CE::element_type)in_data(i) <= in_edges(counts_sz))
    {
      for (j = counts_sz - 1; j >= 0; j--)
      {
        if(in_edges(j) <= (typename CE::element_type)in_data(i))
        {
          counts(j)++;
          if (nargout == 2)
          {
            binidx(i) = j + 1;
          }
          break;
        }
      }
    }
  }
  
  if (nargout == 2)
  {
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    return octave_value(counts);
  }
  
}


template <class CD>
static octave_value_list hco_CountEvenSpacing(CD in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  Matrix binidx;
  octave_idx_type i, j;
  
  double spacing = in_edges(1) - in_edges(0);
  
  #if HCO_DEBUG
    printf("ran hco_CountEvenSpacing\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    binidx = Matrix(in_data.dims());
    binidx.fill(0);
  }
  
  for(i = 0; i < data_sz; i++)
  {
    if((double)in_data(i) == in_edges(counts_sz))
    {
      counts(counts_sz-1)++;
      if (nargout == 2)
      {
        binidx(i) = counts_sz;
      }
    }
    else
    {
      
      j = (octave_idx_type)(((double)in_data(i) - in_edges(0)) / spacing);
      
      // floating points are fickle, so check the output
      if (0 <= j && j <= counts_sz-1)
      {
        
        if ((double)in_data(i) < in_edges(j))
        {
          if(j == 0)
          {
            // data was actually below lowest threshold
            continue;
          }
          else
          {
            j--;
          }
        }
        else if(in_edges(j+1) <= (double)in_data(i))
        {
          if(j == counts_sz-1)
          {
            // data was actually above highest threshold
            continue;
          }
          else
          {
            j++;
          }
        }
        
        counts(j)++;
        
        // should be optimized out
        if (nargout == 2)
        {
          binidx(i) = j+1;
        }
      }
    }
  }
  
  if (nargout == 2)
  {
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    return octave_value(counts);
  }
  
}

template <class CD>
static octave_value_list hco_CountEvenSpacing(CD in_data, FloatMatrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  Matrix binidx;
  octave_idx_type i, j;
  
  float spacing = in_edges(1) - in_edges(0);
    
  #if HCO_DEBUG
    printf("ran hco_CountEvenSpacing\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    binidx = Matrix(in_data.dims());
    binidx.fill(0);
  }
  
  for(i = 0; i < data_sz; i++)
  {
    if((float)in_data(i) == in_edges(counts_sz))
    {
      counts(counts_sz-1)++;
      if (nargout == 2)
      {
        binidx(i) = counts_sz;
      }
    }
    else
    {
      
      j = (octave_idx_type)(((float)in_data(i) - in_edges(0)) / spacing);
      
      // floating points are fickle, so check the output
      if (0 <= j && j <= counts_sz-1)
      {
        
        if ((float)in_data(i) < in_edges(j))
        {
          if(j == 0)
          {
            // data was actually below lowest threshold
            continue;
          }
          else
          {
            j--;
          }
        }
        else if(in_edges(j+1) <= (float)in_data(i))
        {
          if(j == counts_sz-1)
          {
            // data was actually above highest threshold
            continue;
          }
          else
          {
            j++;
          }
        }
        
        counts(j)++;
        
        // should be optimized out
        if (nargout == 2)
        {
          binidx(i) = j+1;
        }
      }
    }
  }
  
  if (nargout == 2)
  {
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    return octave_value(counts);
  }
  
}

template <class CD, class CE>
static octave_value_list hco_CountEvenSpacing(CD in_data, CE in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  Matrix binidx;
  octave_idx_type i, j;
  
  typename CE::element_type spacing = in_edges(1) - in_edges(0);
    
  #if HCO_DEBUG
    printf("ran hco_CountEvenSpacing\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    binidx = Matrix(in_data.dims());
    binidx.fill(0);
  }
  
  for(i = 0; i < data_sz; i++)
  {
    if((typename CE::element_type)in_data(i) == in_edges(counts_sz))
    {
      counts(counts_sz-1)++;
      if (nargout == 2)
      {
        binidx(i) = counts_sz;
      }
    }
    else
    {
      
      j = (octave_idx_type)(((typename CE::element_type)in_data(i) - in_edges(0)) / spacing);
      
      // floating points are fickle, so check the output
      if (0 <= j && j <= counts_sz-1)
      {          
        counts(j)++;
        if (nargout == 2)
        {
          binidx(i) = j+1;
        }
      }
    }
  }
  
  if (nargout == 2)
  {
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    return octave_value(counts);
  }
  
}

template <class CD, class CE>
static octave_value_list hco_CountSorted(CD in_data, CE in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout, sortmode sort_mode)
{
  Matrix binidx;
  octave_idx_type i, j;
  
  #if HCO_DEBUG
    printf("ran hco_CountSorted\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  
  if (nargout == 2)
  {
    binidx = Matrix(in_data.dims());
    binidx.fill(0);
  }
  
  if(sort_mode == ASCENDING)
  {
    for (i = data_sz - 1; i >= 0; i--)
    {
      if ((typename CE::element_type)in_data(i) <= in_edges(counts_sz))
      {
        break;
      }
    }
    
    for (j = counts_sz - 1; j >= 0 && i >= 0; j--)
    {
      for(; i >= 0 && in_edges(j) <= (typename CE::element_type)in_data(i); i--)
      {
        counts(j)++;
        if(nargout == 2)
        {
          binidx(i) = j+1;
        }
      }
    }
    
  }
  else
  {
    for (i = 0; i < data_sz; i++)
    {
      if ((typename CE::element_type)in_data(i) <= in_edges(counts_sz))
      {
        break;
      }
    }
    
    for (j = counts_sz - 1; j >= 0 && i < data_sz; j--)
    {
      for(; i < data_sz && in_edges(j) <= (typename CE::element_type)in_data(i); i++)
      {
        counts(j)++;
        if(nargout == 2)
        {
          binidx(i) = j+1;
        }
      }
    }
  }
  
  if (nargout == 2)
  {
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    return octave_value(counts);
  }
  
}

template <class CD, class CE>
static octave_value_list hco_RunBinningFunc(CD in_data, CE in_edges, sortmode sort_mode, int nargout)
{
  octave_idx_type data_sz = in_data.numel();
  octave_idx_type counts_sz = in_edges.numel()-1;
  
  if (sort_mode || (sort_mode = in_data.issorted()))
  {
    return hco_CountSorted(in_data, in_edges, data_sz, counts_sz, nargout, sort_mode);
  }
  else
  {
    if(hco_CheckEvenSpacing(in_edges, counts_sz))
    {
      return hco_CountEvenSpacing(in_data, in_edges, data_sz, counts_sz, nargout);
    }
    else
    {
      return hco_CountNormal(in_data, in_edges, data_sz, counts_sz, nargout);
    }
  }
}

template <class CD>
static octave_value_list hco_ParseEdgesArg(CD in_data, octave_value in_edges_arg, sortmode sort_mode, int nargout)
{
  if(in_edges_arg.is_double_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.matrix_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.isfloat())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.float_matrix_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_int8_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.int8_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_int16_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.int16_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_int32_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.int32_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_int64_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.int64_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_uint8_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.uint8_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_uint16_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.uint16_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_uint32_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.uint32_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.is_uint64_type())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.uint64_array_value(), sort_mode, nargout);
  }
  else if(in_edges_arg.islogical())
  {
    return hco_RunBinningFunc(in_data, in_edges_arg.int8_array_value(), sort_mode, nargout);
  }
  else
  {
    print_usage ();
    return octave_value();
  }
}


static octave_value_list hco_ParseDataArg(octave_value in_data_arg, octave_value in_edges_arg, sortmode sort_mode, int nargout)
{
  if(in_data_arg.is_double_type())
  {
    return hco_ParseEdgesArg(in_data_arg.matrix_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.isfloat())
  {
    return hco_ParseEdgesArg(in_data_arg.float_matrix_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_int8_type())
  {
    return hco_ParseEdgesArg(in_data_arg.int8_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_int16_type())
  {
    return hco_ParseEdgesArg(in_data_arg.int16_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_int32_type())
  {
    return hco_ParseEdgesArg(in_data_arg.int32_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_int64_type())
  {
    return hco_ParseEdgesArg(in_data_arg.int64_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_uint8_type())
  {
    return hco_ParseEdgesArg(in_data_arg.uint8_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_uint16_type())
  {
    return hco_ParseEdgesArg(in_data_arg.uint16_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_uint32_type())
  {
    return hco_ParseEdgesArg(in_data_arg.uint32_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.is_uint64_type())
  {
    return hco_ParseEdgesArg(in_data_arg.uint64_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else if(in_data_arg.islogical())
  {
    return hco_ParseEdgesArg(in_data_arg.int8_array_value(), in_edges_arg, sort_mode, nargout);
  }
  else
  {
    print_usage ();
    return octave_value();
  }
}


DEFUN_DLD (histcountsoct, args, nargout, "internal function for binning in histcounts")
{
  if(args.length() != 3)
  {
    print_usage ();
    return octave_value(0);
  }
  
  return hco_ParseDataArg(args(0), args(1), (sortmode)args(2).bool_value(), nargout);
  
}
