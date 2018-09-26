#include <octave/oct.h>

DEFUN_DLD (histcountsoct, args, nargout, "internal function for binning in histcounts")
{
    int nargin = args.length();
    octave_value in_data  = args(0);
    octave_value in_edges = args(1);
  
    octave_idx_type i, j;
    octave_idx_type data_sz   = in_data.numel();
    octave_idx_type counts_sz = in_edges.numel()-1;

    octave_value_list ret;

    Matrix counts (1, in_edges.numel()-1);
    counts.fill(0);
    if (nargout == 2)
    {
        Matrix binidx (in_data.dims());
        binidx.fill(0);
        for(i = 0; i < data_sz; i++)
        {
            if(in_data(i) > in_edges(counts_sz))
            {
                continue;
            }
            for(j = counts_sz - 1; j >= 0; j--)
            {
                if(in_edges(j) <= in_data(i))
                {
                    counts(j)++;
                    binidx(i) = j + 1;
                    break;
                }
            }
        }
        ret(1) = octave_value(counts);
        ret(2) = octave_value(binidx);
        return ret;
    }
    else
    {
        for(i = 0; i < data_sz; i++)
        {
            if(in_data(i) > in_edges(counts_sz))
            {
                continue;
            }
            for(j = counts_sz - 1; j >= 0; j--)
            {
                if(in_edges(j) <= in_data(i))
                {
                    counts(j)++;
                  break;
                }
            }
        }
        ret(1) = octave_value(counts);
        return ret;
    }
}
