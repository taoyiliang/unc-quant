#include "Cells.h"

Cell::Cell(std::vector<double> pos_, Material * mtl_, 
          int cid_, int rid_,//linear ID
          double dx_, double dy_,
          int regx_, int regy_, //x,y index of the region
          int ix_, int iy_) //x, y index of cell within region
          : pos(pos_),
            cid(cid_),
            rid(rid_),
            regx(regx_),
            regy(regy_),
            dx(dx_),
            dy(dy_),
            ix(ix_),
            iy(iy_),
            mtl(mtl_){};

//Accessors

double Cell::getTop()   {return pos[1]+0.5*dy;}
double Cell::getBottom(){return pos[1]-0.5*dy;}
double Cell::getLeft()  {return pos[0]+0.5*dx;}
double Cell::getRight() {return pos[0]-0.5*dx;}

//this is the 1g version
/*vector<int> Cell::getIJG()
{
  vector<int> ijg;
  ijg.push_back(ix);
  ijg.push_back(iy);
  return ijg;
}*/

/* //this is the 2g version TODO
vector<int> Cell::getIJG()
{
  vector<int> ijg;
  ijg.push_back(ix);
  ijg.push_back(iy);
  return ijg;
}*/
