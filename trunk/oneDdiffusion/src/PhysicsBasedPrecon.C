#include "PhysicsBasedPrecon.h"

////
// problem specific includes
#include "TwoDProblem.h"

////
// constractor
PhysicsBasedPrecon::PhysicsBasedPrecon(TwoDProblem &problem)
  :_problem(problem)
{}

////
// destractor.
PhysicsBasedPrecon::~PhysicsBasedPrecon()
{}

////
// initialization function
bool PhysicsBasedPrecon::init()
{
  std::cout << "PhysicsBasedPrecon::init()"<<std::endl;
  
  return true;
}

////
// reinitializing preconditioner.
bool PhysicsBasedPrecon::reinit()
{
  std::cout << "PhysicsBasedPrecon::reinit()"<<std::endl;
  return true;
}

int PhysicsBasedPrecon::SetUseTranspose(bool UseTranspose)
{
   /*
    * Disable this option
    */
  return false;
}

int
PhysicsBasedPrecon::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
//  Y=X;
  std::cout << "PhysicsBasedPrecon::apply()"<<std::endl;
  return 0; 
}

////
// getting corresponding problem of interest
TwoDProblem & PhysicsBasedPrecon::getProblem()
{
  return _problem;
}

////
// getting reference to PBP operator
Epetra_Operator & PhysicsBasedPrecon::getPBPOperator()
{
  return *_pbp_operator;
}

double
PhysicsBasedPrecon::NormInf() const
{
   /*
    * Not implemented: Throw an error!
    */
  std::cout << "ERROR: int Preconditioner::NormInf() - "
            << "method is NOT implemented!!  " << std::endl;
  throw "Preconditioner Error";
  return 0.0;
}

const char*
PhysicsBasedPrecon::Label() const
{
  return "No Label is specified in NoPreconditioner...";
  
}

bool
PhysicsBasedPrecon::UseTranspose() const
{
  return false;
  
   /*
    * Not implemented: Throw an error!
    */
  std::cout << "ERROR: Preconditioner::UseTranspose() - "
            << "method is NOT implemented!!  " << std::endl;
  throw "Preconditioner Error";
  return false;
}

bool
PhysicsBasedPrecon::HasNormInf() const
{
   /*
    * NormInf is not implemented
    */
  return false;
}
