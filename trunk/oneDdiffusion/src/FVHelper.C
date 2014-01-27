#include <iostream>

#include "Epetra_Vector.h"
#include "FVHelper.h"
////
// get the neiboring cells in terms of global id, also returns i and j index corresponds to global_id.

namespace FVHelper
{
  void get_neighbors(const int global_id, int &i, int &j, int nx, int ny, int &left, int &right, int &lower, int &upper)
  {
    ////
    // this is "zero" index based..
    // global_id = j*nx+i   (0<=i<=nx, 0<=j<=ny)
    
    i = global_id%nx;
    j = (global_id-i)/nx;

    if(i >0)
      left = global_id-1;
    else
      left = -1;
    if(i<nx-1)
      right = global_id+1;
    else
      right = -1;
    if(j >0 )
      lower = global_id-nx;
    else
      lower = -1;
    if(j < ny-1)
      upper = global_id+nx;
    else
      upper = -1;
    
    return;
  }

  void get_neighbors(const int global_id,int &ieq, int &i, int &j, int nx, int ny, int &left, int &right, int &lower, int &upper)
  {
    ////
    // this is "zero" index based..
    // global_id = j*nx+i   (0<=i<=nx, 0<=j<=ny)
    int nxny = nx*ny;
    int ij = global_id%nxny;

    i = ij%nx;
    ieq = global_id/nxny;
    j = (global_id-i-ieq*nxny)/nx;

    if(i >0)
      left = global_id-1;
    else
      left = -1;
    if(i<nx-1)
      right = global_id+1;
    else
      right = -1;
    if(j >0 )
      lower = global_id-nx;
    else
      lower = -1;
    if(j < ny-1)
      upper = global_id+nx;
    else
      upper = -1;
    
    return;
  }

  void wynn_epsilon(std::vector<double> &sk, int kl, double &sa, int np)
  {
    std::vector<std::vector<double> > skj;
    skj.resize(np);
    for(int i=0;i<np;++i)
      skj[i].resize(np,0);
////
// Aiken
    sa=sk[kl];
    
    int saz=sa;
    double td=1;
    for(int k=0;k<=kl-2;++k)
    {
      skj[k][0]=sk[k];
      double tn=sk[k+2]-sk[k+1];
      td=sk[k+2]-2*sk[k+1]+sk[k];
      
      skj[k][2]=sk[k+2]-(tn*tn)/(td+1.0e-20);
      sa=skj[k][2];
    }
    
    skj[kl][0]=sk[kl];

    for(int j=4;j<=kl;j+=2)
    {
      int jl=j;
      for(int k=0;k<=kl-j;++k)
      {
        double t1=(skj[k+1][j-2]-skj[k][j-2]);
        double t2=(skj[k+1][j-2]-skj[k+2][j-2]);
        double t3=(skj[k+1][j-2]-skj[k+2][j-4]);

        if(fabs(t1*t2*t3)==0)
          return;
        
        skj[k][j]=skj[k+1][j-2]-1.0/(1.0/t1+1.0/t2-1.0/t3);
        double z1=1.0/t1+1.0/t2-1.0/t3;
        sa=skj[k][j];
      }
      
    }
    return;
  }

  void pause()
  {
    std::cin.ignore(1,'\n');
    return;
  }

}
