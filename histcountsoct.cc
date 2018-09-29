#include <octave/oct.h>
#include <limits>

static int hco_CheckEvenSpacing(Matrix in_data, Matrix in_edges, octave_idx_type counts_sz);
static octave_value_list hco_CountEvenSpacing(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout);
static octave_value_list hco_CountNormal(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout);
static octave_value_list hco_CountSorted(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout, sortmode sort_mode);

#define HCO_DEBUG 0

DEFUN_DLD (histcountsoct, args, nargout, "internal function for binning in histcounts")
{
  int nargin = args.length();
  Matrix   in_data   = args(0).matrix_value();
  Matrix   in_edges  = args(1).matrix_value();
  sortmode sort_mode = (sortmode)args(2).double_value();
  
  octave_idx_type i;
  octave_idx_type data_sz   = in_data.numel();
  octave_idx_type counts_sz = in_edges.numel()-1;
  
  if (sort_mode || (sort_mode = in_data.issorted()))
  {
    return hco_CountSorted(in_data, in_edges, data_sz, counts_sz, nargout, sort_mode);
  }
  else
  {
    if(hco_CheckEvenSpacing(in_data, in_edges, counts_sz))
    {
      return hco_CountEvenSpacing(in_data, in_edges, data_sz, counts_sz, nargout);
    }
    else
    {
      return hco_CountNormal(in_data, in_edges, data_sz, counts_sz, nargout);
    }
  }
}

static int hco_CheckEvenSpacing(Matrix in_data, Matrix in_edges, octave_idx_type counts_sz)
{
  octave_idx_type i;
  double spacing = in_edges(1) - in_edges(0);
  for (i = 2; i < counts_sz+1; i++)
  {
    if(abs(spacing - (in_edges(i) - in_edges(i-1))) > std::numeric_limits<double>::epsilon())
    {
      return false;
    }
  }
  return true;
}


static octave_value_list hco_CountEvenSpacing(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  
  octave_idx_type i, j;
  
  double spacing = in_edges(1) - in_edges(0);
    
  #if HCO_DEBUG
    printf("ran hco_CountEvenSpacing\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    Matrix binidx (in_data.dims());
    binidx.fill(0);
    for(i = 0; i < data_sz; i++)
    {
      if(in_data(i) == in_edges(counts_sz))
      {
        counts(counts_sz-1)++;
        binidx(i) = counts_sz;
      }
      else
      {
        
        if(in_edges(0) <= in_data(i) && in_edges < in_data(counts_sz))
        {
          j = (octave_idx_type)((in_data(i) - in_edges(0)) / spacing);
          
          // floating points are fickle, so check the output
          
          if(j == 0)
          {
            if(in_edges(j+1) <= in_data(i))
            {
              j++;
            }
          }
          else if(j == counts_sz-1)
          {
            if(in_data(i) < in_edges(j))
            {
              j--;
            }
          }
          if (0 <= j && j < counts_sz)
          {
            if(in_data(i) < in_edges(j))
            {
              
            }
            else if(in_edges(j+1) <= 
            counts(j)++;
            binidx(i) = j+1;
          }
        }
      }
    }
    octave_value_list ret (2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    for(i = 0; i < data_sz; i++)
    {
      // floating points are fickle, so check the output
      if(in_data(i) == in_edges(counts_sz))
      {
        counts(counts_sz-1)++;
      }
      else
      {
        j = (octave_idx_type)((in_data(i) - in_edges(0)) / spacing);
        
        // floating points are fickle, so check the output
        if (0 <= j && j < counts_sz)
        {
          counts(j)++;
        }
      }
    }
    return octave_value(counts);
  }
}


static octave_value_list hco_CountNormal(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout)
{
  
  octave_idx_type i, j;
  
  #if HCO_DEBUG
    printf("ran hco_CountNormal\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  if (nargout == 2)
  {
    Matrix binidx (in_data.dims());
    binidx.fill(0);
    for(i = 0; i < data_sz; i++)
    {
      if(in_data(i) <= in_edges(counts_sz))
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
    octave_value_list ret(2);
    ret(0) = octave_value(counts);
    ret(1) = octave_value(binidx);
    return ret;
  }
  else
  {
    for(i = 0; i < data_sz; i++)
    {
      if(in_data(i) <= in_edges(counts_sz))
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


static octave_value_list hco_CountSorted(Matrix in_data, Matrix in_edges, octave_idx_type data_sz, octave_idx_type counts_sz, int nargout, sortmode sort_mode)
{
  
  octave_idx_type i, j;
  
  #if HCO_DEBUG
    printf("ran hco_CountSorted\n");
  #endif
  
  Matrix counts (1, counts_sz);
  counts.fill(0);
  
  if(sort_mode == ASCENDING)
  {
    for (i = data_sz - 1; i >= 0; i--)
    {
      if (in_data(i) <= in_edges(counts_sz))
      {
        break;
      }
    }
    
    if (nargout == 2)
    {
      Matrix binidx (in_data.dims());
      binidx.fill(0);
      
      for (j = counts_sz - 1; j >= 0 && i >= 0; j--)
      {
        while(i >= 0 && in_edges(j) <= in_data(i))
        {
          counts(j)++;
          binidx(i) = j+1;
          i--;
        }
      }
      
      octave_value_list ret (2);
      ret(0) = octave_value(counts);
      ret(1) = octave_value(binidx);
      return ret;
    }
    else
    {
      for (j = counts_sz - 1; j >= 0 && i >= 0; j--)
      {
        while(i >= 0 && in_edges(j) <= in_data(i))
        {
          counts(j)++;
          i--;
        }
      }
      return octave_value(counts);
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
    
    if (nargout == 2)
    {
      Matrix binidx (in_data.dims());
      binidx.fill(0);
      
      for (j = counts_sz - 1; j >= 0 && i < data_sz; j--)
      {
        while(i < data_sz && in_edges(j) <= in_data(i))
        {
          counts(j)++;
          binidx(i) = j+1;
          i++;
        }
      }
      
      octave_value_list ret (2);
      ret(0) = octave_value(counts);
      ret(1) = octave_value(binidx);
      return ret;
    }
    else
    {
      for (j = counts_sz - 1; j >= 0 && i < data_sz; j--)
      {
        while(i < data_sz && in_edges(j) <= in_data(i))
        {
          counts(j)++;
          i++;
        }
      }
      return octave_value(counts);
    }
  }
}
