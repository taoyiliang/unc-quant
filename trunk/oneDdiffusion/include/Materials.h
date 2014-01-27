#ifndef MATERIALS_H
#define MATERIALS_H

#include <string>

class Material
{
  public:
    Material(){}
    virtual ~Material(){}
    //vectors for groups
    int                 label;
    std::string         name;
    std::vector<double> siga;
    std::vector<double> nsigf;
    std::vector<double> sigtr;
    std::vector<double> D;
    std::vector<double> chi;

    Material(int                 label_,
             std::string         name_,
             std::vector<double> D_,
             std::vector<double> siga_,
             std::vector<double> nsigf_,
             std::vector<double> sigtr_,
             std::vector<double> chi_):
      label(label_),
      name(name_),
      siga(siga_),
      nsigf(nsigf_),
      sigtr(sigtr_),
      D(D_),
      chi(chi_){};
};


#endif
