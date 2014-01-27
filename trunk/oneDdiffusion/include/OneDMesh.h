#ifndef TWODMESH_H
#define TWODMESH_H

#include <vector>

class TwoDMesh
{
public:
  TwoDMesh();

  virtual ~TwoDMesh()
  {
  }

  void init(int nXregs_, 
            int nYregs_, 
            std::vector<double> &dx,
            std::vector<double> &dy,
            std::vector<int>    &ncells_per_Xreg_,
            std::vector<int>    &ncells_per_Yreg_,
            );
  
  
  int nXregs;// total number of regions in x.(different xs)
  int nYregs;// total number of regions in y.(different xs)
  int ncells;//total number of cells in system
  
  std::vector<int> ncells_per_Xreg;//#of cells per X region,
  std::vector<int> ncells_per_Yreg;//#of cells per Y region,
  std::vector<int> reg_id;// region id for each cell
  
  std::vector<double> x_face; // cell-edge location, ncells+1
  std::vector<double> y_face; // cell-edge location, ncells+1
  std::vector<double> x_center;// cell-center location ncells
  std::vector<double> y_center;// cell-center location ncells
  
};




#endif
