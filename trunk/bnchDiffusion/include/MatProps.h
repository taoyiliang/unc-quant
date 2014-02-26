#ifndef MATPROPS_H
#define MATPROPS_H

#include <vector>

class MatProps
{
public:
  MatProps();

  virtual ~MatProps()
  {
  }


  void init(std::vector<double> &txs_, std::vector<double> &sxs_, std::vector<double> &fxs_, std::vector<double> &nu_);

  std::vector<double> txs;
  std::vector<double> axs;
  std::vector<double> sxs;
  std::vector<double> fxs;
  std::vector<double> nu;
  std::vector<double> nusigf;
  std::vector<double> dfn;
  
};



#endif
