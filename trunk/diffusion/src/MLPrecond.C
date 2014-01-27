#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Ifpack_Utils.h"

#include "TwoDProblem.h"
#include "PhysicsBasedPrecon.h"
#include "MLPrecond.h"

/**
* MLPreconditioner functionality
*/

MLPreconditioner::MLPreconditioner(std::string label,
                                   TwoDProblem &problem)
  :PhysicsBasedPrecon(problem),
   _label(label),
   _problem(problem)
{
//  std::cout << "MLPreconditioner constractor.."<< _label<< std::endl;
  
   /*
    *  Specify Aztec as the smoother.  Here we use subdomain ILU for the preconditioner.
    */

//		Ifpack_PrintSparsity(*M.get(), "matrix.ps", 1);
//		exit(0);

// Set up using Aztec's ILU(0) as the smoother for ML: ie. one level Schwarz with incomplete factorizations
// on each proc, no fill in, no overlap
//
  // set defaults for classic smoothed aggregation with heavy smoothers
    // (of domain decomposition type, i.e. one-level Schwarz with incomplete
    // factorizations on each subdomain/process)
    // We need to define the solvers on each subdomain (== processor).
    // Here we use an incomplete LU factorization, with no fill-in
    // and no overlap. To that aim, we use Aztec's preconditioning function.
    // Aztec requires two more vectors. Note: the following options and params
    // will be used ONLY for the smoother, and will NOT affect the Aztec solver
    // NOTE: to use exact solvers change to AZ_lu (requires AztecOO configured
    // with option--enable-aztecoo-azlu), of use IFPACK smoothers 
    // (requires Trilinos to be built with options--enable-ifpack --enable-amesos)
// 

  // set defaults for classic smoothed aggregation with heavy smoothers
  // (of domain decomposition type, i.e. one-level Schwarz with incomplete
  // factorizations on each subdomain/process)
  // We need to define the solvers on each subdomain (== processor).
  //p Here we use an incomplete LU factorization, with no fill-in
  // and no overlap. To that aim, we use Aztec's preconditioning function.
  // Aztec requires two more vectors. Note: the following options and params
  // will be used ONLY for the smoother, and will NOT affect the Aztec solver
  // NOTE: to use exact solvers change to AZ_lu (requires AztecOO configured
  // with option--enable-aztecoo-azlu), of use IFPACK smoothers 
  // (requires Trilinos to be built with options--enable-ifpack --enable-amesos)

   AZ_defaults(options,params);
   options[AZ_precond] = AZ_dom_decomp;
   options[AZ_subdomain_solve] = AZ_ilu;
   options[AZ_graph_fill] = 0;
   options[AZ_overlap] = 0;

  // SetDefaults() will call AZ_defaults(options,params), and will also set the
  // preconditioner as `AZ_dom_decomp'. 
  // NOTE THAT THE `options' AND `params' VECTORS ARE NOT COPIED into
  // the list, only the pointers is stored, so do not delete options 
  // and params before the end of the linear system solution!
  // Alternatively, you can also call SetDefaults() without passing 
  // `options' and `params.' This way, the code will allocate a int 
  // and a double vector, that must be freed by the user.
  // `DD' means to set default values for domain decomposition
  // preconditioners


//   ML_Epetra::SetDefaults("SA", precParams, options, params);
   ML_Epetra::SetDefaults("SA", precParams, options, params);

  // Overwrite some parameters. Please refer to the user's guide
  // for more information
  // Some parameters are reported here to better explain the process
  // even if they are as defaults. 
  // NOTE: To use `METIS' as aggregation scheme, you need to configure
  // ML with the option --with-ml_metis. Otherwise, the code will
  // creates aggregates containing all the local nodes (that is,
  // the dimension of the coarse problem will be equal to the
  // number of processors)


//   Use Aztec for a smoother
   precParams.set("smoother: type", "Aztec");
//  precParams.set("smoother: type", "Chebyshev");
   precParams.set("smoother: sweeps", 5);

  // Put 64 nodes on each aggregate. This number can be too small
  // for large problems. In this case, either augment this value, or increase
  // the number of levels. Also, use only presmoothing, and KLU as
  // coarse solver (KLU is enabled by default with Amesos)


//   precParams.set("aggregation: type", "METIS");
//   precParams.set("aggregation: nodes per aggregate", 128);
//   precParams.set("aggregation: type", "Uncoupled-MIS");

// Direct solve on the coarse mesh
   precParams.set("coarse: type", "Amesos-KLU");

///   precParams.set("prec type", "full-MGV");
   precParams.set("max levels", 1);
   
  // output level
//  precParams.set("ML output", 20);
  /*
  // max number of level
  precParams.set("max levels", 5);
  // set finest level to 0
  precParams.set("increasing or decreasing","increasing");

  // use uncoupled scheme to create the aggregate
  precParams.set("aggregation: type", "Uncoupled");

  ////
  // Smoothing option
  // Chebyshev smoother.
  precParams.set("smoother: type", "Chebyshev");
  // number of sweeps
  precParams.set("smoother: sweeps",1);
  // use both pre and post smoothing.
  precParams.set("smoother: pre or post", "both");

#ifdef HAVE_ML_AMESOS
  //solve with serial direct solver KLU
  precParams.set("coarse: type", "Amesos-KLU");
#else
  //if no amsos, use jacobi...
  precParams.set("coarse: type", "Jacobi");
#endif  
   
  precParams.set("PDE equations", 1);
  */
//   precParams.set("ML print initial list", 0);

//  precParams.set("smoother: type","symmetric Gauss-Seidel");

  // use both pre and post smoothing
  precParams.set("smoother: pre or post", "both");

  ////
  // allocate preconditioning matrix.
  if( _M !=NULL )
    delete _M;
  

//  if( _M == NULL)
    _M = new Epetra_CrsMatrix(Copy, _problem.getMap(), 5);
}


MLPreconditioner::~MLPreconditioner()
{
  //delete all the preconditioning matrix
  delete _Prec;

  delete _M;
  
}

void
//MLPreconditioner::initialize( const Teuchos::RCP<Epetra_CrsMatrix>& M, Teuchos::ParameterList* parameters )
MLPreconditioner::initialize()
{
 std::cout << "MLPreconditioner::initialize()"<<std::endl;
 throw;
 
 _Prec->ComputePreconditioner();
}

void
MLPreconditioner::setup(const Epetra_Vector *x)
{
 std::cout << "MLPreconditioner::setup()"<<std::endl;

  ////
  // setup matrix.
 _problem.createMatrix(x,_M);
  
 if( _Prec != NULL)
   delete _Prec;
    
 _Prec = new ML_Epetra::MultiLevelPreconditioner( *_M, precParams, true);
  
}

int
MLPreconditioner::ApplyInverse( const Epetra_MultiVector& r, Epetra_MultiVector& z ) const
{
   /*
    * pass along the request to ML.
    */
//  std::cout << r;
//  //  std::cout << z;

  //z=r;
  //return 0;

  int ok=_Prec->ApplyInverse(r,z);
  
  return 0;

}

const Epetra_Comm&
MLPreconditioner::Comm() const
{
  return _problem.getComm();
}

const Epetra_Map& 
MLPreconditioner::OperatorDomainMap () const
{
  ////
  // domain map appear to be same as standard map.
  return _problem.getMap();
  
}

const Epetra_Map& 
MLPreconditioner::OperatorRangeMap () const
{
  ////
  // range map appear to be same as standard map.
  return _problem.getMap();
  
}

std::string
MLPreconditioner::getLabel()
{
  return _label;
  
}
