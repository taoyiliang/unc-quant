#ifndef _PHYSICSBASEDPRECON_H
#define _PHYSICSBASEDPRECON_H

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"

////
// problem specific includes

class TwoDProblem;

// Physics-based preconditioner Class
class PhysicsBasedPrecon : public Epetra_Operator
{ 

public:

  ////
  // Constructor
  PhysicsBasedPrecon(TwoDProblem &problem);

  ////
  // Destractor
  virtual ~PhysicsBasedPrecon();

  ////
  // initialization
  bool init();

  ////
  // reinitialization
  bool reinit();

  virtual int SetUseTranspose(bool UseTranspose);
  
  virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
  
  virtual double NormInf() const;
  
  virtual const char* Label() const;
  
  virtual bool UseTranspose() const;
  
  virtual bool HasNormInf() const;

  ////
  // function to get reference to the problem of interest
  TwoDProblem & getProblem();
  

  ////
  // getting PBP operator
  Epetra_Operator & getPBPOperator();
  
private:
  
  TwoDProblem &_problem;
  Epetra_Operator *_pbp_operator;
  
  
};


#endif //_PHYSICSBASEDPRECON_H
