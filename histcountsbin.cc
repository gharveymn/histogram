#include <limits>
#include <octave/oct.h>

static bool
hco_CheckEvenSpacing (NDArray in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  double spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz + 1; i++)
    {
      if (fabs (spacing - (in_edges(i) - in_edges(i - 1)))
          > std::numeric_limits<double>::epsilon ())
        return false;
    }
  return true;
}

static bool
hco_CheckEvenSpacing (FloatNDArray in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  float spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz + 1; i++)
    {
      if (fabsf (spacing - (in_edges(i) - in_edges(i - 1)))
          > std::numeric_limits<float>::epsilon ())
        return false;
    }
  return true;
}

template <class C>
static bool
hco_CheckEvenSpacing (C in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  typename C::element_type spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz + 1; i++)
    {
      if (spacing != (in_edges(i) - in_edges(i - 1)))
        return false;
    }
  return true;
}

template <class CD, class CE>
static octave_value_list
hco_CountNormal (CD in_data, CE in_edges, octave_idx_type data_sz,
                 octave_idx_type counts_sz, int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

  Matrix counts (1, counts_sz);
  counts.fill (0);
  if (nargout == 2)
    {
      binidx = NDArray (in_data.dims ());
      binidx.fill (0);
    }

  for (i = 0; i < data_sz; i++)
    {
      if (in_data(i) <= in_edges(counts_sz))
        {
          for (j = counts_sz - 1; j >= 0; j--)
            {
              if (in_edges(j) <= in_data(i))
                {
                  counts(j)++;
                  if (nargout == 2) // should be optimized out
                    binidx(i) = j + 1;
                  break;
                }
            }
        }
    }

  if (nargout == 2)
    {
      octave_value_list ret (2);
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    return octave_value (counts);
}

template <class CD>
static octave_value_list
hco_CountEvenSpacing (CD in_data, NDArray in_edges, octave_idx_type data_sz,
                      octave_idx_type counts_sz, int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

  double spacing = in_edges(1) - in_edges(0);

  Matrix counts (1, counts_sz);
  counts.fill (0);
  if (nargout == 2)
    {
      binidx = NDArray (in_data.dims ());
      binidx.fill (0);
    }

  for (i = 0; i < data_sz; i++)
    {
      if (in_data(i) == in_edges(counts_sz))
        {
          counts (counts_sz - 1)++;
          if (nargout == 2)
            binidx(i) = counts_sz;
        }
      else if (in_edges(0) <= in_data(i) && in_data(i) < in_edges(counts_sz))
        {
          j = (octave_idx_type)(((double)in_data(i) - in_edges(0)) / spacing);
          // floating points are fickle, so check the output
          if (j == -1)
            j = 0;
          else if (j == counts_sz)
            j = counts_sz-1;
          else if (in_data(i) < in_edges(j))
            j--;
          else if (in_edges(j + 1) <= in_data(i))
            j++;
          counts(j)++;
          if (nargout == 2) // should be optimized out
            binidx(i) = j + 1;
        }
    }

  if (nargout == 2)
    {
      octave_value_list ret (2);
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    return octave_value (counts);
}

template <class CD>
static octave_value_list
hco_CountEvenSpacing (CD in_data, FloatNDArray in_edges,
                      octave_idx_type data_sz, octave_idx_type counts_sz,
                      int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

  float spacing = in_edges(1) - in_edges(0);

  Matrix counts (1, counts_sz);
  counts.fill (0);
  if (nargout == 2)
    {
      binidx = NDArray (in_data.dims ());
      binidx.fill (0);
    }

  for (i = 0; i < data_sz; i++)
    {
      if (in_data(i) == in_edges(counts_sz))
        {
          counts(counts_sz - 1)++;
          if (nargout == 2)
            binidx(i) = counts_sz;
        }
      else if (in_edges(0) <= in_data(i) && in_data(i) < in_edges(counts_sz))
        {
          j = (octave_idx_type)(((float)in_data(i) - in_edges(0)) / spacing);
          // floating points are fickle, so check the output
          if (j == -1)
            j = 0;
          else if (j == counts_sz)
            j = counts_sz-1;
          else if (in_data(i) < in_edges(j))
            j--;
          else if (in_edges(j + 1) <= in_data(i))
            j++;
          counts(j)++;
          if (nargout == 2) // should be optimized out
            binidx(i) = j + 1;
        }
    }

  if (nargout == 2)
    {
      octave_value_list ret (2);
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    return octave_value (counts);
}

template <class CD, class CE>
static octave_value_list
hco_CountEvenSpacing (CD in_data, CE in_edges, octave_idx_type data_sz,
                      octave_idx_type counts_sz, int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

  typename CE::element_type::val_type spacing = in_edges(1) - in_edges(0);

  Matrix counts (1, counts_sz);
  counts.fill (0);
  if (nargout == 2)
    {
      binidx = NDArray (in_data.dims ());
      binidx.fill (0);
    }

  for (i = 0; i < data_sz; i++)
    {
      if (in_data(i) == in_edges(counts_sz))
        {
          counts (counts_sz - 1)++;
          if (nargout == 2)
            binidx(i) = counts_sz;
        }
      else if (in_edges(0) <= in_data(i) && in_data(i) < in_edges(counts_sz))
        {
          j = (octave_idx_type) (
              (typename CE::element_type::val_type)(
              (typename CE::element_type)in_data(i) - in_edges(0))
              / spacing);
          // printf("j: %li, (cast)in_data: %f, in_edges(0): %f, sub: %f, div: %f, spacing: %f\n", j, (double)((typename CE::element_type)in_data(i)), (double)(in_edges(0)), (double)((typename CE::element_type)in_data(i) - in_edges(0)), (double)((typename CE::element_type::val_type)((typename CE::element_type)in_data(i) - in_edges(0)) / spacing), (double)spacing);
          counts(j)++;
          if (nargout == 2)
            binidx(i) = j + 1;
        }
    }

  if (nargout == 2)
    {
      octave_value_list ret (2);
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    return octave_value (counts);
}

template <class CD, class CE>
static octave_value_list
hco_CountSorted (CD in_data, CE in_edges, octave_idx_type data_sz,
                 octave_idx_type counts_sz, int nargout, sortmode sort_mode)
{
  NDArray binidx;
  octave_idx_type i, j;

  Matrix counts (1, counts_sz);
  counts.fill (0);

  if (nargout == 2)
    {
      binidx = NDArray (in_data.dims ());
      binidx.fill (0);
    }

  if (sort_mode == ASCENDING)
    {
      for (i = data_sz - 1; i >= 0; i--)
        {
          if (in_data(i) <= in_edges(counts_sz))
            break;
        }
      for (j = counts_sz - 1; j >= 0 && i >= 0; j--)
        {
          for (; i >= 0 && in_edges(j) <= in_data(i); i--)
            {
              counts(j)++;
              if (nargout == 2)
                binidx(i) = j + 1;
            }
        }
    }
  else
    {
      for (i = 0; i < data_sz; i++)
        {
          if (in_data(i) <= in_edges(counts_sz))
            break;
        }
      for (j = counts_sz - 1; j >= 0 && i < data_sz; j--)
        {
          for (; i < data_sz && in_edges(j) <= in_data(i); i++)
            {
              counts(j)++;
              if (nargout == 2)
                binidx(i) = j + 1;
            }
        }
    }

  if (nargout == 2)
    {
      octave_value_list ret (2);
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    return octave_value (counts);
}

template <class CD, class CE>
static octave_value_list
hco_RunBinningFunc (CD in_data, CE in_edges, sortmode sort_mode, int nargout)
{
  octave_idx_type data_sz = in_data.numel ();
  octave_idx_type counts_sz = in_edges.numel () - 1;

  if (sort_mode || (sort_mode = in_data.issorted ()))
    return hco_CountSorted (in_data, in_edges, data_sz, counts_sz, nargout,
                            sort_mode);
  else
    {
      if (hco_CheckEvenSpacing (in_edges, counts_sz))
        return hco_CountEvenSpacing (in_data, in_edges, data_sz, counts_sz,
                                     nargout);
      else
        return hco_CountNormal (in_data, in_edges, data_sz, counts_sz,
                                nargout);
    }
}

template <class CD>
static octave_value_list
hco_ParseEdgesArg (CD in_data, octave_value in_edges_arg, sortmode sort_mode,
                   int nargout)
{
  if (in_edges_arg.is_double_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.isfloat ())
    return hco_RunBinningFunc (in_data, in_edges_arg.float_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_int8_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.int8_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_int16_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.int16_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_int32_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.int32_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_int64_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.int64_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_uint8_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.uint8_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_uint16_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.uint16_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_uint32_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.uint32_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.is_uint64_type ())
    return hco_RunBinningFunc (in_data, in_edges_arg.uint64_array_value (),
                               sort_mode, nargout);
  else if (in_edges_arg.islogical ())
    return hco_RunBinningFunc (in_data, in_edges_arg.int8_array_value (),
                               sort_mode, nargout);
    
  print_usage ();
  return octave_value ();
}

static octave_value_list
hco_ParseDataArg (octave_value in_data_arg, octave_value in_edges_arg,
                  sortmode sort_mode, int nargout)
{
  if (in_data_arg.is_double_type ())
    return hco_ParseEdgesArg (in_data_arg.array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.isfloat ())
    return hco_ParseEdgesArg (in_data_arg.float_array_value (), in_edges_arg, 
                              sort_mode, nargout);
  else if (in_data_arg.is_int8_type ())
    return hco_ParseEdgesArg (in_data_arg.int8_array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.is_int16_type ())
    return hco_ParseEdgesArg (in_data_arg.int16_array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.is_int32_type ())
    return hco_ParseEdgesArg (in_data_arg.int32_array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.is_int64_type ())
    return hco_ParseEdgesArg (in_data_arg.int64_array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.is_uint8_type ())
    return hco_ParseEdgesArg (in_data_arg.uint8_array_value (), in_edges_arg,
                              sort_mode, nargout);
  else if (in_data_arg.is_uint16_type ())
    return hco_ParseEdgesArg (in_data_arg.uint16_array_value (),
                              in_edges_arg, sort_mode, nargout);
  else if (in_data_arg.is_uint32_type ())
    return hco_ParseEdgesArg (in_data_arg.uint32_array_value (),
                              in_edges_arg, sort_mode, nargout);
  else if (in_data_arg.is_uint64_type ())
    return hco_ParseEdgesArg (in_data_arg.uint64_array_value (),
                              in_edges_arg, sort_mode, nargout);
  else if (in_data_arg.islogical ())
    return hco_ParseEdgesArg (in_data_arg.int8_array_value (), in_edges_arg,
                              sort_mode, nargout);
  print_usage ();
  return octave_value ();
}

DEFUN_DLD (histcountsbin, args, nargout,
           "internal function for binning in histcounts")
{
  if (args.length () != 3)
    {
      print_usage ();
      return octave_value (0);
    }

  return hco_ParseDataArg (args(0), args(1),
                           (sortmode)args(2).bool_value (), nargout);
}
