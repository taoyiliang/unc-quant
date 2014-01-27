
/** \class Preconditioner
 *
 * Preconditioner Class.  Subclass of
 * Epetra_Operator.
 */

#include <string>

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "Ifpack.h"

#ifndef included_preconditioner
#define included_preconditioner

class PhysicsPrecond : public Epetra_Operator
{
public:
   /**
    * Constructor. 
    */
  PhysicsPrecond(string label, const Epetra_Map &map); 

   /**
    * Destructor.
    */
   virtual ~PhysicsPrecond();

//   virtual void initialize( const Teuchos::RCP<Epetra_CrsMatrix>& M, Teuchos::ParameterList* parameters ) = 0;
//   virtual void initialize() ;
//   virtual void initialize() = 0;

   virtual int SetUseTranspose(bool UseTranspose);

   virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

   /**
    * Apply the preconditioner.  This is currently just a solve with ML.
    */

   virtual int ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

   virtual double NormInf() const;

   virtual const char* Label() const;

   virtual bool UseTranspose() const;

   virtual bool HasNormInf() const;

   virtual const Epetra_Comm& Comm() const; // Must be supplied by children

   virtual const Epetra_Map& OperatorDomainMap() const ; // Must be supplied by children

   virtual const Epetra_Map& OperatorRangeMap() const;  // Must be supplied by children

protected:

   string d_label;
   const Epetra_Map & _my_map;
   
   
};

#endif
