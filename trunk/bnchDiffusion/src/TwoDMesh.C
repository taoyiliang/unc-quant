#include "TwoDMesh.h"
#include <iostream>

TwoDMesh::TwoDMesh()
{
}

int TwoDMesh::cid(int x, int y)
{
  return x + y*nX;
}

int TwoDMesh::gid(int x, int y, int g)
{
  return x + y*nX + g*nX*nY;
}

int TwoDMesh::getG(int gid)
{
  return ijg(gid)[2];
}

int TwoDMesh::getI(int gid)
{
  return ijg(gid)[0];
}

int TwoDMesh::getJ(int gid)
{
  return ijg(gid)[1];
}

std::vector<int> TwoDMesh::ij(int cid)
{
  //this ignores group and returns cell x,y
  std::vector<int> coords;
  int x = cid % nX;
  int y = (cid - x)/nX;
  coords.push_back(x);
  coords.push_back(y);
  return coords;
}

// Uncomment this when you go to 2D_2g
std::vector<int> TwoDMesh::ijg(int gid)
{
  //this is the 2D_2g implementation
  std::vector<int> coords;
  int nXY = nX*nY;
  int temp = gid % nXY;

  int g = (gid - temp)/nXY;
  int x = temp % nX;
  int y = (temp - x)/nX;

  coords.push_back(x);
  coords.push_back(y);
  coords.push_back(g);
  return coords;
}

void TwoDMesh::init(int nXregs_,
                    int nYregs_,
                    double &dx, 
                    double &dy, 
                    int &ncells_per_Xreg_,
                    int &ncells_per_Yreg_,
                    std::vector<int> &matid)
{
  nXregs = nXregs_;
  nYregs = nYregs_;
  nCells_per_Xreg = ncells_per_Xreg_;
  nCells_per_Yreg = ncells_per_Yreg_;
  nX = nCells_per_Xreg * nXregs;
  nY = nCells_per_Yreg * nYregs;
/*  std::cout<<"In TwoDMesh::init"<<std::endl;
  std::cout<<"n_regs: "<<nXregs_<<" "<<nYregs_<<std::endl;
  std::cout<<"dx/dy : "<<dx<<" "<<dy<<std::endl;
  std::cout<<"n/reg : "<<ncells_per_Xreg_<<" "<<ncells_per_Yreg_<<std::endl;*/
  std::cout<<"matid : "<<std::endl;
  for (size_t i=0;i<matid.size();++i)
  {
    if (i%11==0 and i!=0)
      std::cout<<std::endl;
    std::cout<<matid[i]<<" ";
  }
  std::cout<<std::endl;
  ncells = 0;
  double x(0.0),y(0.0);
  std::vector<double> loc(2,0.0);

  Material * curMat;
  Cell * newCell;

  int xi,yi,xregi,global_id,region_id;
  int yregi = -1;
  for (int iYreg=0;iYreg<nYregs;++iYreg)
  {
    yregi+=1;
    yi = yregi*nCells_per_Yreg - 1;
    loc[1] = 0.5*dy + yregi*nCells_per_Yreg;
    for (int iY=0;iY<nCells_per_Yreg;++iY)
    {
      yi+=1;
      loc[1]=(0.5+yi)*dy;
      xregi=-1;
      for (int iXreg=0;iXreg<nXregs;++iXreg)
      {
        xregi+=1;
        region_id=nXregs*iYreg + iXreg;
        xi = xregi*nCells_per_Xreg - 1;
        for (int iX=0;iX<nCells_per_Xreg;++iX)
        {
          xi+=1;
          loc[0]=(0.5+xi)*dx;
          global_id=cid(xi,yi);

          curMat = &materials[matid[region_id]];
          newCell = new Cell(loc,&*curMat,global_id,region_id,dx,dy,xregi,yregi,xi,yi);
          //cells.push_back(&(new Cell(loc,&*curMat,global_id,dx,dy,
          //                         xregi,yregi,xi,yi)));
          cells.push_back(*newCell);
          ncells+=1;
        }//iX
      }//iXreg
    }//iY
  }//iYreg

  //DEBUG print cell stuff
  /*std::cout<<"Cell gid:"<<std::endl;
  for (int gij=0;gij<nX*nY;++gij)
  {
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    if (gij < 10)
      std::cout<<" ";
    if (gij < 100)
      std::cout<<" ";
    std::cout<<cells[gij].gid<<" ";
  }
  std::cout<<std::endl;*/
  //DEBUG print cell stuff
/*  std::cout<<"Cell ij:"<<std::endl;
  int myi,myj;
  for (int gij=0;gij<nX*nY;++gij)
  {
    myi = cells[gij].ix;
    myj = cells[gij].iy;
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    if (myi < 10)
      std::cout<<" ";
    std::cout<<cells[gij].ix<<",";
    if (myj < 10)
      std::cout<<" ";
    std::cout<<cells[gij].iy<<" ";
  }
  std::cout<<std::endl;
  //DEBUG print cell stuff
  std::cout<<"Cell matID:"<<std::endl;
  for (int gij=0;gij<nX*nY;++gij)
  {
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    std::cout<<cells[gij].mtl->name<<" ";
  }
  std::cout<<std::endl;
*/
  //DEBUG print cell stuff
  /*std::cout<<"Cell D:"<<std::endl;
  for (int gij=0;gij<nX*nY;++gij)
  {
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    std::cout<<cells[gij].mtl->D[1]<<" ";
  }
  std::cout<<std::endl;*/
  //DEBUG print cell stuff
  /*std::cout<<"Cell up:"<<std::endl;
  for (int gij=0;gij<nX*nY;++gij)
  {
    int jp,i,gijp;
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    if (cells[gij].iy==0)
      std::cout<<" -  ";
    else
    {
      i = cells[gij].ix;
      jp = cells[gij].iy-1;
      gijp = gid(i,jp);
    if (gijp < 10)
      std::cout<<" ";
    if (gijp < 100)
      std::cout<<" ";
      std::cout<<cells[gijp].gid<<" ";
    }
  }
  std::cout<<std::endl;*/
  //DEBUG print cell stuff
  /*std::cout<<"Cell down:"<<std::endl;
  for (int gij=0;gij<nX*nY;++gij)
  {
    int jm,i,gijp;
    if (gij%nX==0 and gij!=0)
      std::cout<<std::endl;
    if (cells[gij].iy==nY-1)
      std::cout<<" -  ";
    else
    {
      i = cells[gij].ix;
      jm = cells[gij].iy+1;
      gijp = gid(i,jm);
    if (gijp < 10)
      std::cout<<" ";
    if (gijp < 100)
      std::cout<<" ";
      std::cout<<cells[gijp].gid<<" ";
    }
  }
  std::cout<<std::endl;*/
}//TwoDMesh::init

void TwoDMesh::add_material(int label, 
                            std::string name,
                            std::vector<double> D,
                            std::vector<double> siga, 
                            std::vector<double> nsigf,
                            std::vector<double> sigtr, 
                            std::vector<double> chi)
{
  materials[label]=*(new Material(label,name,D,siga,nsigf,sigtr,chi));
  std::cout<<"added material "<<materials[label].label<<": "
                              <<materials[label].name<<std::endl;
}
