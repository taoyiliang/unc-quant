#include "PhysicsPrecond.h"



PhysicsPrecond::PhysicsPrecond( string label, const Epetra_Map & map)
  :_my_map(map),
  d_label(label)
{
  
  std::cout <<"exsiting physicsprecond constractor"<<std::endl;
}


PhysicsPrecond::~PhysicsPrecond()
{
}

int
PhysicsPrecond::SetUseTranspose(bool UseTranspose)
{
   /*
    * Disable this option
    */
  return false;
}

int
PhysicsPrecond::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   /*
    * Not implemented: Throw an error!
    */
  cout << "ERROR: Preconditioner::Apply() - "
       << "method is NOT implemented!!  " << endl;
  throw "Preconditionerr Error";
  return false;
}

int
PhysicsPrecond::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  std::cout << "PhysicsPrecond::ApplyInverse: " << std::endl;
  
  Y=X;
  std::cout << "returning.."<<std::endl;
  
  return 0;
  
   /*
    * Not implemented: Throw an error!
    */
  cout << "ERROR: Preconditioner::Apply() - "
       << "method is NOT implemented!!  " << endl;
  throw "Preconditionerr Error";
  return false;
}

double
PhysicsPrecond::NormInf() const
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
PhysicsPrecond::Label() const
{
   return d_label.c_str();
}

bool
PhysicsPrecond::UseTranspose() const
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
PhysicsPrecond::HasNormInf() const
{
   /*
    * NormInf is not implemented
    */
  return false;
}


const Epetra_Comm& PhysicsPrecond::Comm() const
{
  return _my_map.Comm();
}


const Epetra_Map& PhysicsPrecond::OperatorDomainMap() const
{
  return _my_map;
}

const Epetra_Map& PhysicsPrecond::OperatorRangeMap() const
{
  return _my_map;
}
