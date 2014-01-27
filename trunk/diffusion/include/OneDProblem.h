#ifndef TWODSNPROBLEM_H
#define TWODSNPROBLEM_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

#include "NOX.H"
#include "NOX_Epetra.H"

#include "TwoDMesh.h"
// Forward Declarations
class Epetra_Import;
class Epetra_Operator;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

class MatProps;
class TwoDMaps;

// Flag to tell the evaluate routine what objects to fill
enum FillType {FILL_RES, FILL_JAC, FILL_PREC, FILL_FD_RES,FILL_MF_RES,FILL_MF_JAC,FILL_USER}; 
enum BCs {VACUUM, REFLECTIVE};

class TwoDProblem
{
  
public:
  TwoDProblem(TwoDMesh &mesh, MatProps &matpro, Epetra_Comm &comm );


  ~TwoDProblem()
  {
  }


  bool evaluate(FillType f, const Epetra_Vector *solnVector, 
		Epetra_Vector *rhsVector);


  void createMatrix(const Epetra_Vector *x, Epetra_CrsMatrix *_M);
  
  void evaluateResidual(const Epetra_Vector* soln, Epetra_Vector *tmp_rhs);
  
  double computeCurrent(int icell, const Epetra_Vector * soln);
  
    
  bool setupPrecOperator(const Epetra_Vector *x, Epetra_Operator *Prec, Teuchos::ParameterList* precParams);

  double computeKeff(const Epetra_Vector *soln);

  double computeFissionSourceStrength(double keff, const Epetra_Vector *soln);
  
  void copyOldSolution(Epetra_Vector &soln);

  void setBCs(int left, int right);
  
  ////-------
  // accessor.
  Teuchos::RCP<Epetra_Vector> getSolution()
  {
    return _new_soln;
  }

  Teuchos::RCP<Epetra_Vector> getOldSolution()
  {
    return _old_soln;
  }

  Epetra_Comm &getComm()
  {
    return _comm;
  }
  
  Epetra_Map &getMap()
  {
    return *_map;
  }

  TwoDMesh &getMesh()
  {
    return _mesh;
  }

  bool setKeff(double keff)
  {
    _keff = keff;
  }

  Teuchos::RCP<Epetra_CrsMatrix> GetMatrix()
  {
    return _Jac;
  }
  
private:
  Epetra_Map *_map;
  Epetra_Map *_face_map;
  
  Epetra_Comm &_comm;

  Teuchos::RCP<Epetra_Vector> _new_soln;
  Teuchos::RCP<Epetra_Vector> _old_soln;

  TwoDMesh &_mesh;
  MatProps &_matprop;
  
  bool _flag;

  BCs _left_bc;
  BCs _right_bc;
  
  
  double _keff;
  double _keff_old;


  Teuchos::RCP<Epetra_CrsMatrix> _Jac;
  
};

#endif
