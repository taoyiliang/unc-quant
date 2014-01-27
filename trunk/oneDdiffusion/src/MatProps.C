#include <iostream>
#include "MatProps.h"

MatProps::MatProps()
{
}


void MatProps::init(std::vector<double> &txs_, std::vector<double> &sxs_, std::vector<double> &fxs_, std::vector<double> &nu_)
{
  for(int ireg=0;ireg<txs_.size();++ireg)
  {
    txs.push_back(txs_[ireg]);
    sxs.push_back(sxs_[ireg]);
    fxs.push_back(fxs_[ireg]);
    nu.push_back(nu_[ireg]);

    
    double _nusigf = nu[ireg]*fxs[ireg];
    double _axs    = txs[ireg]-sxs[ireg];
    double _dfn    = 1/(3*txs[ireg]);

    nusigf.push_back(_nusigf);
    axs.push_back(_axs);
    dfn.push_back(_dfn);
  }
  
  return;
}

