#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Ifpack_Utils.h"

#include "TwoDProblem.h"
#include "PhysicsBasedPrecon.h"
#include "NoPreconditioner.h"

/**
* Example of concrete implementation of physics-based preconditioner,
* this class does not have acutral preconditioning step, just operate with identity matrix.
*/
  

NoPreconditioner2::NoPreconditioner2( TwoDProblem &problem)
  :PhysicsBasedPrecon(problem),
   _problem(problem)
{
  std::cout << "constractor of NoPreconditioner.."<<std::endl;
  
}


NoPreconditioner2::~NoPreconditioner2()
{}

int
NoPreconditioner2::ApplyInverse( const Epetra_MultiVector& r, Epetra_MultiVector& z ) const
{
   /*
    * just use identity operator.
    */
//  std::cout << "NoPreconditioner::ApplyInverse"<< r.NumVectors()<<std::endl;
//  std::cout << r;
  
  z=r;
  double result;
  
  z.Norm2(&result);
//  std::cout << "z-norm: "<<result<<std::endl;
  return 0;
}

double
NoPreconditioner2::NormInf() const
{
   /*
    * Not implemented: Throw an error!
    */
  cout << "ERROR: int Preconditioner::NormInf() - "
       << "method is NOT implemented!!  " << endl;
  throw "Preconditioner Error";
  return 0.0;
}

const char*
NoPreconditioner2::Label() const
{
  return "No Label is specified in NoPreconditioner...";
  
}

bool
NoPreconditioner2::UseTranspose() const
{
  return false;
  
   /*
    * Not implemented: Throw an error!
    */
  cout << "ERROR: Preconditioner::UseTranspose() - "
       << "method is NOT implemented!!  " << endl;
  throw "Preconditioner Error";
  return false;
}

bool
NoPreconditioner2::HasNormInf() const
{
   /*
    * NormInf is not implemented
    */
  return false;
}

const Epetra_Comm&
NoPreconditioner2::Comm() const
{
  return _problem.getComm();
}

const Epetra_Map& 
NoPreconditioner2::OperatorDomainMap () const
{
  return _problem.getMap();
}

const Epetra_Map& 
NoPreconditioner2::OperatorRangeMap () const
{
  return _problem.getMap();
}
