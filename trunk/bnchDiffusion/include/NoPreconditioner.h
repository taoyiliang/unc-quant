
/** \class Preconditioner
 *
 * Preconditioner Class.  Subclass of
 * Epetra_Operator.
 */


#ifndef NOPRECONDITIONER_H
#define NOPRECONDITIONER_H

////
// forward declarations
class PhysicsBasedPrecon;


class NoPreconditioner2 : public PhysicsBasedPrecon
{
public:
   /**
    * Constructor. 
    */
  NoPreconditioner2( TwoDProblem &problem ); 

   /**
    * Destructor.
    */
   virtual ~NoPreconditioner2();

   /**
    * Apply the preconditioner.  This is currently just a solve with identity matrix.
    */
   
   virtual int ApplyInverse( const Epetra_MultiVector& r, Epetra_MultiVector& z ) const;

   virtual double NormInf() const;

   virtual const char* Label() const;

   virtual bool UseTranspose() const;

   virtual bool HasNormInf() const;

   virtual const Epetra_Comm& Comm() const;

   virtual const Epetra_Map& OperatorDomainMap() const;

   virtual const Epetra_Map& OperatorRangeMap() const;

private:

   TwoDProblem &_problem;
   
   Teuchos::ParameterList precParams;
   
};

#endif
