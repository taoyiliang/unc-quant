#include "NOX_Common.H"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#include "TwoDProblem.h"

#include "MatProps.h"
#include "MLPrecond.h"
#include "NeutronicsParam.h"

TwoDProblem::TwoDProblem(TwoDMesh &mesh, Epetra_Comm &comm ) :
  _mesh(mesh),
  _comm(comm),
  _keff(1),
  _left_bc(VACUUM),
  _right_bc(VACUUM),
  _top_bc(VACUUM),
  _bottom_bc(VACUUM),
  _nca_flag(false)
{
  //std::cout<<"ncells: "<<_mesh.ncells<<std::endl;
  _map      = new Epetra_Map(_mesh.ncells*2,0,comm);
  _face_map = new Epetra_Map(_mesh.ncells*2+1,0,comm);

  _new_soln = Teuchos::rcp(new Epetra_Vector(*_map));
  _new_soln->PutScalar(0.01);

  _old_soln = Teuchos::rcp(new Epetra_Vector(*_new_soln));

  _Jac = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *_map, 5));
  
}

//-----------------------------------------------------------------------------------------
bool TwoDProblem::evaluate(FillType f, const Epetra_Vector* soln, Epetra_Vector* tmp_rhs)
{

  std::cout << "evaluate"<<std::endl;
  
  _flag = f;

  if ( _nca_flag )
    _keff_nl = computeNonlinearKeff(*soln);
  
  ////
  // LO residual function evaluation.
  evaluateResidual(soln, tmp_rhs);

  
  return true;

}

//-------------------------------------------------------------------------------------
void TwoDProblem::evaluateResidual(const Epetra_Vector* soln, Epetra_Vector *tmp_rhs)
{
  //std::cout<<"In TwoDProblem::evaluateResidual:"<<std::endl;
  //std::cout<<"i,S,abs,DX,DY"<<endl;
  for (int g=0;g<2;++g)
  {
    for(int icell=0;icell< _mesh.ncells;++icell)
    {
      std::vector<int> ij=_mesh.ij(icell);
      int i=ij[0];
      int j=ij[1];

      std::vector<int> gid_ctr;
      gid_ctr.push_back(_mesh.gid(i,j,0));
      gid_ctr.push_back(_mesh.gid(i,j,1));
      double source(0);
      double removal = _mesh.cells[icell].mtl->siga[g]*(*soln)[gid_ctr[g]];
      double scatt = _mesh.cells[icell].mtl->sigtr[0]*(*soln)[gid_ctr[0]];
 
      if (g==0)
      {
        removal += scatt;
        if (_nca_flag)
        {
          source =  _mesh.cells[icell].mtl->nsigf[0]*(*soln)[gid_ctr[0]]/_keff_nl;
          source += _mesh.cells[icell].mtl->nsigf[1]*(*soln)[gid_ctr[1]]/_keff_nl;
        }
        else
        {
          source =  _mesh.cells[icell].mtl->nsigf[0]*(*_old_soln)[gid_ctr[0]]/_keff;
          source += _mesh.cells[icell].mtl->nsigf[1]*(*_old_soln)[gid_ctr[1]]/_keff;
        }
      }
      else if (g==1)
      {
        source = scatt;
      }
      else
      {
        std::cout<<"This is a problem!  g not in (0,1): "<<g<<std::endl;
        std::cin.get();
      }

      double currentX = computeCurrentX(icell,g,soln);
      double currentY = computeCurrentY(icell,g,soln);
      double current = currentX+currentY;
      //std::cout<<"diff: "<<g<<" "<<i<<" "<<j<<" "<<current<<std::endl;
      (*tmp_rhs)[gid_ctr[g]] =  currentX + currentY + removal -source;
      //std::cout<<"end   g,i,j: "<<g<<" "<<i<<" "<<j<<std::endl;
      //std::cout<<"g,i,j|cX,cY,r,s: "<<g<<" "<<i<<" "<<j<<" | ";
      //std::cout<<currentX<<" "<<currentY<<" "<<removal<<" "<<source<<std::endl;

      std::string fileName1("curX.out");
      std::string fileName2("curY.out");
      std::string fileName3("remv.out");
      std::string fileName4("srcs.out");
      //AppendOutVals(i,j,g,currentX,fileName1);
      //AppendOutVals(i,j,g,currentY,fileName2);
      //AppendOutVals(i,j,g,removal,fileName3);
      //AppendOutVals(i,j,g,source,fileName4);
    }//icell
  }

  return;
}

void AppendOutVals(int i, int j, int g,double val,std::string fileName)
{
  std::ofstream out;
  out.precision(10);
  out<<std::scientific;
  out.open(fileName.c_str(), std::ios_base::app);

  out<<i<<" "<<j<<" "<<g<<" "<<val<<std::endl;
}

//----------------------------------------------------------
double TwoDProblem::computeKeff(const Epetra_Vector *soln)
{
  ////
  // assuming (1/k)*nsigf[0]*phi is normalized to one.
  double fission_src = 0;
  for(int icell =0;icell<_mesh.ncells;++icell)
  {
    int i,j;
    std::vector<int> ij=_mesh.ij(icell);
    i=ij[0];
    j=ij[1];
    double dx = _mesh.cells[icell].dx;
    double dy = _mesh.cells[icell].dy;
    for (int g=0;g<2;++g)
    {
      int curGID=_mesh.gid(i,j,g);
      fission_src += _mesh.cells[icell].mtl->nsigf[g]*(*soln)[curGID]*dx*dy;
    }
  }
  return fission_src;
  
}

//----------------------------------------------------------
double TwoDProblem::computeNonlinearKeff(const Epetra_Vector &soln)
{
  double fission_src = 0;
  for ( int icell = 0; icell<_mesh.ncells;++icell)
  {
    int i,j;
    std::vector<int> ij=_mesh.ij(icell);
    i=ij[0];
    j=ij[1];
    double dx = _mesh.cells[icell].dx;
    double dy = _mesh.cells[icell].dy;
    for (int g=0;g<2;++g)
    {
      int curGID=_mesh.gid(i,j,g);
      fission_src += _mesh.cells[icell].mtl->nsigf[g]*soln[curGID]*dx*dy;
    }
  }
  return fission_src;
}

//-------------------------------------------------------
double TwoDProblem::computeFissionSourceStrength(double keff, const Epetra_Vector *soln)
{
  ////
  // assuming (1/k)*nsigf[0]*phi is normalized to one.
  double fission_src = 0;
  for(int icell =0;icell<_mesh.ncells;++icell)
  {
    int i,j;
    std::vector<int> ij=_mesh.ij(icell);
    i=ij[0];
    j=ij[1];
    double dx = _mesh.cells[icell].dx;
    double dy = _mesh.cells[icell].dy;
    for (int g=0;g<2;++g)
    {
      int curGID=_mesh.gid(i,j,g);
      fission_src += _mesh.cells[icell].mtl->nsigf[g]*(*soln)[curGID]*dx*dy;
    }
  }
  return fission_src/keff;
  
}

//--------------------------------------------------------------------------
double TwoDProblem::computeCurrentX(int icell, int g, const Epetra_Vector * soln)
{
  double right_phi, left_phi;
  
  double dx = _mesh.cells[icell].dx;
  //double dy = _mesh.cells[icell].dy;
  //double beta = dx/dy;
  double dx_p,dx_m;

  Cell * cellC;
  Cell * cellR;
  Cell * cellL;

  cellC = &(_mesh.cells[icell]);
  int i = cellC->ix;
  int j = cellC->iy;
  int gid_right,gid_left;
  int cid_right,cid_left;
  int gid_ctr = _mesh.gid(i,j,g);

  //int ireg = _mesh.reg_id[icell];
  //int reg_r,reg_l;

  double dfn_term = 0;

  if( i == 0 )
  {
    gid_right = _mesh.gid(i+1,j,g);
    cid_right = _mesh.cid(i+1,j);
    cellR = &(_mesh.cells[cid_right]);

    double dx_r  = cellR->dx;
    //double dfn_r = (cellR->mtl->D[g]*dx_r + cellC->mtl->D[g]*dx)/(dx_r+dx);
    double dfn_r = 0.5*(cellR->mtl->D[g] + cellC->mtl->D[g]);
    //dx_r = (dx_r + dx)*0.5;
    double ddphidx_r = dfn_r*((*soln)[gid_ctr]-(*soln)[gid_right])/dx/dx;
    
    double dx_l  = dx;
    double dfn_l = cellC->mtl->D[g];
    double left_phi = -(*soln)[gid_ctr]*(0.25-dfn_l/dx)/(0.25+dfn_l/dx);//vacuum
    double ddphidx_l = dfn_l*((*soln)[gid_ctr]-left_phi)/dx/dx;
      
    double current = 0.0;

    if( _left_bc == VACUUM)
      current = ddphidx_l;
    
    current = ddphidx_r;
    
    
    return current;
    
  }
  else if (i == _mesh.nX-1)
  {
    int cid_left = _mesh.cid(i-1,j);
    int gid_left = _mesh.gid(i-1,j,g);
    cellL = &(_mesh.cells[cid_left]);

    //int reg_l = _mesh.reg_id[icell-1];
    double dx_l = cellL->dx;
    double dfn_l = (cellL->mtl->D[g]*dx_l + cellC->mtl->D[g]*dx)/(dx_l+dx);
    dx_l = (dx_l + dx)*0.5;
    double ddphidx_l = dfn_l*((*soln)[gid_ctr]-(*soln)[gid_left])/dx_l;
      
    double dx_r = dx;
    //int reg_r = _mesh.reg_id[icell];
    double dfn_r = cellC->mtl->D[g];
    double right_phi = -(*soln)[gid_ctr]*(0.25-dfn_r/dx)/(0.25+dfn_r/dx);//vacuum
    double ddphidx_r = dfn_r*((*soln)[gid_ctr]-right_phi)/dx_r;
      
    double current = 0.0;

    if( _right_bc == VACUUM )
      current = ddphidx_r;
    
    current += ddphidx_l;
    
    return current/dx;
  }
  else
  {
  int cid_left = _mesh.cid(i-1,j);
  int gid_left = _mesh.gid(i-1,j,g);
  cellL = &(_mesh.cells[cid_left]);
  int cid_right = _mesh.cid(i+1,j);
  int gid_right = _mesh.gid(i+1,j,g);
  cellR = &(_mesh.cells[cid_right]);

    //int reg_l = _mesh.reg_id[icell-1];
    double dx_l  = cellL->dx;
    //double dfn_l = (cellL->mtl->D[g]*dx_l + cellC->mtl->D[g]*dx)/(dx_l+dx);
    double dfn_l = 0.5*(cellL->mtl->D[g] + cellC->mtl->D[g]);
    //dx_l = (dx_l + dx)*0.5;
    double ddphidx_l = dfn_l*((*soln)[gid_ctr]-(*soln)[gid_left])/dx/dx;
      

    //int reg_r = _mesh.reg_id[icell+1];
    double dx_r = cellR->dx;
    //double dfn_r = (cellR->mtl->D[g]*dx_r + cellC->mtl->D[g]*dx)/(dx_r+dx);
    double dfn_r = 0.5*(cellR->mtl->D[g] + cellC->mtl->D[g]);
    //dx_r = (dx_r + dx)*0.5;
    double ddphidx_r = dfn_r*((*soln)[gid_ctr]-(*soln)[gid_right])/dx/dx;
    
    double current = ddphidx_l + ddphidx_r;
      
    return current;
  }
  
}

//--------------------------------------------------------------------------
double TwoDProblem::computeCurrentY(int icell, int g, const Epetra_Vector * soln)
{
  double top_phi, bot_phi;
  
  double dy = _mesh.cells[icell].dy;
  //double dx = _mesh.cells[icell].dx;
  double dy_p,dy_m;

  Cell * cellC;
  Cell * cellT;
  Cell * cellB;

  cellC = &(_mesh.cells[icell]);
  int i = cellC->ix;
  int j = cellC->iy;
  int gid_ctr = _mesh.gid(i,j,g);
  int gid_top,gid_bottom;
  int cid_top,cid_bottom;

  //int ireg = _mesh.reg_id[icell];
  //int reg_r,reg_l;

  double dfn_term = 0;

  if( j == 0 )
  {
    cid_top = _mesh.cid(i,j+1);
    gid_top = _mesh.gid(i,j+1,g);
    cellT = &(_mesh.cells[cid_top]);

    //int reg_r = _mesh.reg_id[icell+1];
    double dy_t  = cellT->dy;
    double dfn_t = (cellT->mtl->D[g]*dy_t + cellC->mtl->D[g]*dy)/(dy_t+dy);
    dy_t = (dy_t + dy)*0.5;
    double ddphidy_t = dfn_t*((*soln)[gid_ctr]-(*soln)[gid_top])/dy_t;

    
    //int reg_l = _mesh.reg_id[icell];
    double dy_b  = dy;
    double dfn_b = cellC->mtl->D[g];//_matprop.dfn[ireg];
    double bottom_phi = -(*soln)[gid_ctr]*(0.25-dfn_b/dy)/(0.25+dfn_b/dy);//vacuum
    double ddphidy_b = dfn_b*((*soln)[gid_ctr]-bottom_phi)/dy_b;
      
    double current = 0.0;

    if( _bottom_bc == VACUUM)
      current = ddphidy_b;
    
    current += ddphidy_t;
    
    
    return current/dy;
    
  }
  else if (j == _mesh.nY-1)
  {
    int cid_bottom = _mesh.cid(i,j-1);
    int gid_bottom = _mesh.gid(i,j-1,g);
    cellB = &(_mesh.cells[cid_bottom]);

    //int reg_l = _mesh.reg_id[icell-1];
    double dy_b = cellB->dy;
    double dfn_b = (cellB->mtl->D[g]*dy_b + cellC->mtl->D[g]*dy)/(dy_b+dy);
    dy_b = (dy_b + dy)*0.5;
    double ddphidy_b = dfn_b*((*soln)[gid_ctr]-(*soln)[gid_bottom])/dy_b;
      
    double dy_t = dy;
    //int reg_r = _mesh.reg_id[icell];
    double dfn_t = cellC->mtl->D[g];

    double top_phi = -(*soln)[gid_ctr]*(0.25-dfn_t/dy)/(0.25+dfn_t/dy);//vacuum
    double ddphidy_t = dfn_t*((*soln)[gid_ctr]-top_phi)/dy_t;
      
    double current = 0.0;

    if( _top_bc == VACUUM )
      current = ddphidy_t;
    
    current += ddphidy_b;
    
    return current/dy;
  }
  else
  {
    int cid_bottom = _mesh.cid(i,j-1);
    int gid_bottom = _mesh.gid(i,j-1,g);
    cellB = &(_mesh.cells[cid_bottom]);
    int cid_top = _mesh.cid(i,j+1);
    int gid_top = _mesh.gid(i,j+1,g);
    cellT = &(_mesh.cells[cid_top]);
    //int reg_l = _mesh.reg_id[icell-1];
    double dy_b  = cellB->dy;
    double dfn_b = (cellB->mtl->D[g]*dy_b + cellC->mtl->D[g]*dy)/(dy_b+dy);
    dy_b = (dy_b + dy)*0.5;
    double ddphidy_b = dfn_b*((*soln)[gid_ctr]-(*soln)[gid_bottom])/dy_b;
      

    //int reg_r = _mesh.reg_id[icell+1];
    double dy_t = cellT->dy;
    double dfn_t = (cellT->mtl->D[g]*dy_t + cellC->mtl->D[g]*dy)/(dy_t+dy);
    dy_t = (dy_t + dy)*0.5;
    double ddphidy_t = dfn_t*((*soln)[gid_ctr]-(*soln)[gid_top])/dy_t;
    
    double current = ddphidy_b + ddphidy_t;
      
    return current/dy;
  }
  
}

//----------------------------------------------------------------------------------------------------------------------
bool TwoDProblem::setupPrecOperator(const Epetra_Vector *x, Epetra_Operator *Prec, Teuchos::ParameterList* precParams)
{
  
  MLPreconditioner *MLPrec = dynamic_cast<MLPreconditioner *>(Prec);
  MLPrec->setup(x);

  return true;
}

//----------------------------------------------------------------------------
void TwoDProblem::createMatrix(const Epetra_Vector *x, Epetra_CrsMatrix *_M)
{
  
  bool DEBUG = false;
  Cell * cellC;
  Cell * cellT;
  Cell * cellB;
  Cell * cellL;
  Cell * cellR;

  Epetra_Vector ones(*x);

  ones.PutScalar(1.0);
  //ones.Print(std::cout);

  for (int g=0;g<2;++g)
  {
    for(int icell=0;icell<_mesh.ncells;++icell)
    {
      int i,j;
      i = _mesh.cells[icell].ix;
      j = _mesh.cells[icell].iy;
      std::vector<int> gid_ctr;
      gid_ctr.push_back(_mesh.gid(i,j,0));
      gid_ctr.push_back(_mesh.gid(i,j,1));
      //int gid_ctr=_mesh.gid(i,j,g);

      const int gid_bottom = _mesh.gid( i ,j-1,g);
      const int gid_top =    _mesh.gid( i ,j+1,g);
      const int gid_left =   _mesh.gid(i-1, j ,g);
      const int gid_right =  _mesh.gid(i+1, j ,g);

      const int cid_bottom = _mesh.cid( i ,j-1);
      const int cid_top =    _mesh.cid( i ,j+1);
      const int cid_left =   _mesh.cid(i-1, j );
      const int cid_right =  _mesh.cid(i+1, j );

      cellC = &(_mesh.cells[icell]);
      cellB = &(_mesh.cells[cid_bottom]);
      cellT = &(_mesh.cells[cid_top]);
      cellL = &(_mesh.cells[cid_left]);
      cellR = &(_mesh.cells[cid_right]);

      double left_phi, right_phi;
      double bottom_phi, top_phi;

      double dx = cellC->dx;
      double dy = cellC->dy;
      int    nX = _mesh.nX;
      int    nY = _mesh.nY;

      double dfnL,dfnR;
      double dfnC;
      double dfnB,dfnT;


      //removal, source
      double src(0);
      double abs   = _mesh.cells[icell].mtl->siga[g];
      double scatt = _mesh.cells[icell].mtl->sigtr[0]; 
      if (g==0)
      {
        abs += scatt;
        src += cellC->mtl->nsigf[0]/_keff + cellC->mtl->nsigf[1]/_keff;
      }
      else if (g==1)
      {
        src = scatt; //downscatter is dependent on group 0 flux
      }
      
      std::vector<double> value(5);
      std::vector<int> indices(5);
      //now get the current
      double currentX = computeCurrentX(icell,g,&ones);
      double currentY = computeCurrentX(icell,g,&ones);

      double diag = currentX + currentY + abs;
      bool fail_flag = _M->SumIntoGlobalValues(gid_ctr[g],1,&diag,&gid_ctr[g]);
      if(fail_flag)
        _M->InsertGlobalValues(gid_ctr[g],1,&diag,&gid_ctr[g]);
    }
  }
  _M->FillComplete();
  return;
}
/*
      if( i == 0 )
      {
        dfnR = -cellR->mtl->D[g]/dx/dx;
        dfnL = cellC->mtl->D[g];
        dfnC = -dfnR;

        if (j == 0)//Left Bottom
        {
          if (DEBUG) std::cout<<"Left Bottom"<<std::endl;
          dfnT = -cellT->mtl->D[g]/dy/dy;
          dfnB = cellC->mtl->D[g];
          dfnC += -dfnT;

          if (_left_bc == VACUUM)
          {
            dfnC+=dfnL*(1+(0.25-dfnL/dx)/(0.25+dfnL/dx))/dx/dx;
          }
          if (_bottom_bc == VACUUM)
          {
            dfnC+=dfnB*(1+(0.25-dfnB/dy)/(0.25+dfnB/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
        else if (j == _mesh.nY-1)//Left Top
        {
          if (DEBUG) std::cout<<"Left Top"<<std::endl;
          dfnT = cellC->mtl->D[g];
          dfnB = -cellB->mtl->D[g]/dy/dy;
          dfnC += -dfnB;

          if (_left_bc == VACUUM)
          {
            dfnC+=dfnL*(1+(0.25-dfnL/dx)/(0.25+dfnL/dx))/dx/dx;
          }
          if (_top_bc == VACUUM)
          {
            dfnC+=dfnT*(1+(0.25-dfnT/dy)/(0.25+dfnT/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        }
        else //Left Interior
        {
      if (DEBUG) std::cout<<"Left Interior"<<std::endl;
          dfnT = -cellT->mtl->D[g]/dy/dy;
          dfnB = -cellB->mtl->D[g]/dy/dy;
          dfnC += -(dfnT+dfnB);

          if (_left_bc == VACUUM)
          {
            dfnC+=dfnL*(1+(0.25-dfnL/dx)/(0.25+dfnL/dx))/dx/dx;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
      }
      else if (i == _mesh.nX-1)
      {
        dfnR = cellC->mtl->D[g];
        dfnL = -cellL->mtl->D[g]/dx/dx;
        dfnC = -dfnL;

        if (j == 0)//Right Bottom
        {
      if (DEBUG) std::cout<<"Right Bottom"<<std::endl;
          dfnT = -cellT->mtl->D[g]/dy/dy;
          dfnB = cellC->mtl->D[g];
          dfnC += -dfnT;

          if (_right_bc == VACUUM)
          {
            dfnC+=dfnR*(1+(0.25-dfnR/dx)/(0.25+dfnR/dx))/dx/dx;
          }
          if (_bottom_bc == VACUUM)
          {
            dfnC+=dfnB*(1+(0.25-dfnB/dy)/(0.25+dfnB/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
        else if (j == _mesh.nY-1)//Right Top
        {
      if (DEBUG) std::cout<<"Right Top"<<std::endl;
          dfnT = cellC->mtl->D[g];
          dfnB = -cellB->mtl->D[g]/dy/dy;
          dfnC += -dfnB;

          if (_right_bc == VACUUM)
          {
            dfnC+=dfnR*(1+(0.25-dfnR/dx)/(0.25+dfnR/dx))/dx/dx;
          }
          if (_top_bc == VACUUM)
          {
            dfnC+=dfnT*(1+(0.25-dfnT/dy)/(0.25+dfnT/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        }
        else//Right Interior
        {
      if (DEBUG) std::cout<<"Right Interior"<<std::endl;
          dfnT = -cellT->mtl->D[g]/dy/dy;
          dfnB = -cellB->mtl->D[g]/dy/dy;
          dfnC += -(dfnT+dfnB);

          if (_right_bc == VACUUM)
          {
            dfnC+=dfnR*(1+(0.25-dfnR/dx)/(0.25+dfnR/dx))/dx/dx;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
      }
      else
      {
        dfnR = -cellR->mtl->D[g]/dx;
        dfnL = -cellL->mtl->D[g]/dx;
        dfnC = -(dfnR+dfnL)/dx;
        if (j == 0)//Bottom Interior
        {
      if (DEBUG) std::cout<<"Bottom Interior"<<std::endl;
          dfnT = -cellT->mtl->D[g]/dy/dy;
          dfnB = cellC->mtl->D[g];
          dfnC += -dfnT;

          if (_bottom_bc == VACUUM)
          {
            dfnC+=dfnB*(1+(0.25-dfnB/dy)/(0.25+dfnB/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        
          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
        else if (j == _mesh.nY-1)//Top Interior
        {
      if (DEBUG) std::cout<<"Top Interior"<<std::endl;
          dfnT = cellC->mtl->D[g];
          dfnB = -cellB->mtl->D[g]/dy/dy;
          dfnC += -dfnB;

          if (_top_bc == VACUUM)
          {
            dfnC+=dfnT*(1+(0.25-dfnT/dy)/(0.25+dfnT/dy))/dy/dy;
          }

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        
          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        }
        else//Interior Interior
        {
      if (DEBUG) std::cout<<"Total Interior"<<std::endl;
          dfnT =  -cellT->mtl->D[g]/dy/dy;
          dfnB =  -cellB->mtl->D[g]/dy/dy;
          dfnC += -(dfnT+dfnB);

          bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnC,&gid_ctr);

          //int gid_ctrB = gid_ctr-nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnB,&gid_bottom);

          //int gid_ctrL = gid_ctr-1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnL,&gid_left);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnL,&gid_left);
        
          //int gid_ctrR = gid_ctr+1;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnR,&gid_right);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnR,&gid_right);
        
          //int gid_ctrT = gid_ctr+nX;
          fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&dfnT,&gid_top);
          if (fail_flag)
            fail_flag = _M->InsertGlobalValues(gid_ctr,1,&dfnT,&gid_top);
        }
      }
  */    
      /*double fission_term = 0.0;
      double absorption_term = _matprop.axs[ireg];
      double sum_term = -fission_term + absorption_term;
      */
/*      bool fail_flag = _M->SumIntoGlobalValues(gid_ctr,1,&abs, &gid_ctr);
      
    }//icell
  }//g
  _M->FillComplete();
  
  return;
}*/


//------------------------------------------------------
void TwoDProblem::copyOldSolution(Epetra_Vector &soln)
{
  _old_soln->Update(1.0,soln,0.0);

  return;
}


//
void TwoDProblem::setBCs(int left, int right, int bottom, int top)
{
  if( left == 0 )
    _left_bc = VACUUM;
  else
    _left_bc = REFLECTIVE;

  if( right==0 )
    _right_bc = VACUUM;
  else
    _right_bc = REFLECTIVE;

  if( bottom == 0 )
    _bottom_bc = VACUUM;
  else
    _bottom_bc = REFLECTIVE;

  if( top==0 )
    _top_bc = VACUUM;
  else
    _top_bc = REFLECTIVE;

  return;
  
}
