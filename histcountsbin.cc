#include <limits>
#include <octave/oct.h>

#define HCO_DEBUG 0

static bool
hco_CheckEvenSpacing (NDArray in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  double spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz + 1; i++)
    {
      if (fabs (spacing - (in_edges(i) - in_edges(i - 1)))
          > std::numeric_limits<double>::epsilon ())
        {
          return false;
        }
    }
  return true;
}

static octave_value_list
hco_CountNormal (NDArray in_data, NDArray in_edges, octave_idx_type data_sz,
                 octave_idx_type counts_sz, int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

#if HCO_DEBUG
  printf ("ran hco_CountNormal\n");
#endif

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
      ret(0) = octave_value (counts);
      ret(1) = octave_value (binidx);
      return ret;
    }
  else
    {
      return octave_value (counts);
    }
}

static octave_value_list
hco_CountEvenSpacing (NDArray in_data, NDArray in_edges,
                      octave_idx_type data_sz, octave_idx_type counts_sz,
                      int nargout)
{
  NDArray binidx;
  octave_idx_type i, j;

  float spacing = in_edges(1) - in_edges(0);

#if HCO_DEBUG
  printf ("ran hco_CountEvenSpacing\n");
#endif

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
            {
              binidx(i) = counts_sz;
            }
        }
      else
        {
          j = (octave_idx_type)((in_data(i) - in_edges(0)) / spacing);
          // floating points are fickle, so check the output
          if (0 <= j && j <= counts_sz - 1)
            {
              if (in_data(i) < in_edges(j))
                {
                  if (j == 0)
                    {
                      // data was actually below lowest threshold
                      continue;
                    }
                  else
                    {
                      j--;
                    }
                }
              else if (in_edges(j + 1) <= in_data(i))
                {
                  if (j == counts_sz - 1)
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
              if (nargout == 2) // should be optimized out
                {
                  binidx(i) = j + 1;
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
    {
      return octave_value (counts);
    }
}

static octave_value_list
hco_CountSorted (NDArray in_data, NDArray in_edges, octave_idx_type data_sz,
                 octave_idx_type counts_sz, int nargout, sortmode sort_mode)
{
  NDArray binidx;
  octave_idx_type i, j;

#if HCO_DEBUG
  printf ("ran hco_CountSorted\n");
#endif

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
            {
              break;
            }
        }

      for (j = counts_sz - 1; j >= 0 && i >= 0; j--)
        {
          for (; i >= 0
                 && in_edges(j) <= in_data(i);
               i--)
            {
              counts(j)++;
              if (nargout == 2)
                {
                  binidx(i) = j + 1;
                }
            }
        }
    }
  else
    {
      for (i = 0; i < data_sz; i++)
        {
          if (in_data(i) <= in_edges(counts_sz))
            {
              break;
            }
        }

      for (j = counts_sz - 1; j >= 0 && i < data_sz; j--)
        {
          for (; i < data_sz
                 && in_edges(j) <= in_data(i);
               i++)
            {
              counts(j)++;
              if (nargout == 2)
                {
                  binidx(i) = j + 1;
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
    {
      return octave_value (counts);
    }
}

DEFUN_DLD (histcountsbin, args, nargout,
           "internal function for binning in histcounts")
{
  if (args.length () != 3)
    {
      print_usage ();
      return octave_value (0);
    }
    
  NDArray in_data   = args(0).array_value ();
  NDArray in_edges  = args(1).array_value ();
  sortmode sort_mode = (sortmode)args(2).bool_value ();

  octave_idx_type data_sz = in_data.numel ();
  octave_idx_type counts_sz = in_edges.numel () - 1;

  if (sort_mode || (sort_mode = in_data.issorted ()))
    {
      return hco_CountSorted (in_data, in_edges, data_sz, counts_sz, nargout,
                              sort_mode);
    }
  else
    {
      if (hco_CheckEvenSpacing (in_edges, counts_sz))
        {
          return hco_CountEvenSpacing (in_data, in_edges, data_sz, counts_sz,
                                       nargout);
        }
      else
        {
          return hco_CountNormal (in_data, in_edges, data_sz, counts_sz,
                                  nargout);
        }
    }
}
