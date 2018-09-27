#include <octave/oct.h>

DEFUN_DLD (histcountsoct, args, nargout, "internal function for binning in histcounts")
{
    int nargin = args.length();
    Matrix in_data  = args(0).matrix_value();
    Matrix in_edges = args(1).matrix_value();
  
    octave_idx_type i, j;
    octave_idx_type data_sz   = in_data.numel();
    octave_idx_type counts_sz = in_edges.numel()-1;

    Matrix counts (1, counts_sz);
    counts.fill(0);
  
    // check if bin sizes are uniform
  
    if (nargout == 2)
    {
        Matrix binidx (in_data.dims());
        binidx.fill(0);
        for(i = 0; i < data_sz; i++)
        {
            if(in_edges(counts_sz) <= in_data(i))
            {
              for (j = counts_sz - 1; j >= 0; j--)
              {
                  if(in_edges(j) <= in_data(i))
                  {
                      counts(j)++;
                      binidx(i) = j + 1;
                      break;
                  }
              }
            }
        }
        octave_value_list ret;
        ret(0) = octave_value(counts);
        ret(1) = octave_value(binidx);
        return ret;
    }
    else
    {
        for(i = 0; i < data_sz; i++)
        {
            if(in_edges(counts_sz) <= in_data(i))
            {
              for(j = counts_sz - 1; j >= 0; j--)
              {
                  if(in_edges(j) <= in_data(i))
                  {
                    counts(j)++;
                    break;
                  }
              }
            }
        }
        return octave_value(counts);
    }
}
