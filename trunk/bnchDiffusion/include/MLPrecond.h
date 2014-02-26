////
// trilinos includes
#include "Epetra_CrsMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

////
// local includes
#include "PhysicsBasedPrecon.h"

#ifndef included_MLpreconditioner
#define included_MLpreconditioner

class MLPreconditioner : public PhysicsBasedPrecon
{
public:
   /**
    * Constructor. 
    */
  MLPreconditioner( std::string label,
                    TwoDProblem &problem); 

   /**
    * Destructor.
    */
   virtual ~MLPreconditioner();

   /**
    * Initialize from parameter list.
    */
//   void initialize( const Teuchos::RCP<Epetra_CrsMatrix>& M, Teuchos::ParameterList* parameters );
   void initialize();

   /**
    * setting up matrix.
    */
   void setup(const Epetra_Vector *x);
   
   
   /**
    * Apply the preconditioner.  This is currently just a solve with ML.
    */
   virtual int ApplyInverse( const Epetra_MultiVector& r, Epetra_MultiVector& z ) const;

   virtual const Epetra_Comm& Comm() const;

   virtual const Epetra_Map& OperatorDomainMap() const;

   virtual const Epetra_Map& OperatorRangeMap() const;

   std::string getLabel();
   
   Epetra_CrsMatrix *_M;   
private:

   TwoDProblem &_problem;

   ML_Epetra::MultiLevelPreconditioner* _Prec;   

   Teuchos::ParameterList precParams;

   int options[AZ_OPTIONS_SIZE];
   double params[AZ_PARAMS_SIZE];


   int _neq;
   
   std::string _label;
   
};

#endif
