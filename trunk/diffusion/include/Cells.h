#ifndef CELLS_H
#define CELLS_H

#include <vector>

#include "Materials.h"

class Cell
{
  public:
    Cell(){};
    Cell(std::vector<double> pos_, 
         Material * mtl_,
         int cid_, 
         int rid_, 
         double dx_, 
         double dy_, 
         int regx_,
         int regy_,
         int ix_, 
         int iy_);
    virtual ~Cell(){};

    //accessors
    double getTop();
    double getBottom();
    double getLeft();
    double getRight();
    double getIJG();

    //variables
    std::vector<double> pos;
    Material * mtl;
    int cid;
    int rid;
    int regx;
    int regy;
    double dx;
    double dy;
    int ix;
    int iy;
};


#endif
