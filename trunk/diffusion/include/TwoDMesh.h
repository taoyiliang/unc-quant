#ifndef TWODMESH_H
#define TWODMESH_H

#include <map>
#include <vector>
#include <string>
#include "Cells.h"
#include "Materials.h"

class TwoDMesh
{
public:
  TwoDMesh();

  virtual ~TwoDMesh()
  {
  }

  int cid(int,int);
  int gid(int,int,int);
  int getG(int);
  int getI(int);
  int getJ(int);
  std::vector<int> ij(int);
  std::vector<int> ijg(int);

  void init(int nXregs_, 
            int nYregs_, 
            double &dx,
            double &dy,
            int &ncells_per_Xreg,
            int &ncells_per_Yreg,
            std::vector<int> &matid
            );

  void add_material(int,
                    std::string,
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>);

  int getLowerCell(int);
  int getUpperCell(int);
  int getLeftCell(int);
  int getRightCell(int);
  
  int nCells_per_Xreg;
  int nCells_per_Yreg;
  int nXregs;// total number of regions in x.(different xs)
  int nYregs;// total number of regions in y.(different xs)
  int nX; //total cells in X
  int nY; //total cells in Y
  int ncells;//total number of cells in system
  
  std::vector<Cell> cells;

  std::map<int,Material> materials;

  //std::vector<int> ncells_per_Xreg;//#of cells per X region,
  //std::vector<int> ncells_per_Yreg;//#of cells per Y region,
  //std::vector<int> reg_id;// region id for each cell
  
  //std::map<int,std::vector<int,int> > get_2D;//use linear to get x,y ID
  //std::map<std::vector<int,int>,int > get_1D;//use x,y ID to get linear ID
  //std::vector<double> xy_face; // cell-edge location, ncells+1
  //std::vector<double> xy_center;// cell-center location ncells
  
};




#endif
