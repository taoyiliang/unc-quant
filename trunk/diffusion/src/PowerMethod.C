#include "PowerMethod.h"

#include "NeutronicsParam.h"
#include "FVHelper.h"

PowerMethod::PowerMethod(TwoDProblem &problem, Teuchos::RCP<NOX::Solver::Generic> solver, InputOutput &input_output) 
  :  IterativeScheme(problem,solver,input_output)
{
}

PowerMethod::~PowerMethod()
{
}

void PowerMethod::init()
{
}

void PowerMethod::setup()
{
}

bool PowerMethod::run( )
{
  std::cout << "in PowerMethod::run()"<<std::endl;
  double &keff  = NeutronicsParam::keff;
  int    &ngrps = NeutronicsParam::n_prompt_groups;
  int ncells = _problem.getMesh().ncells;
  
  double keff_tol   = _my_input_file("Execution/keff_tol",1e-5);
  int max_power_its = _my_input_file("Execution/max_power_its",2000);
  int max_s_its     = _my_input_file("Execution/max_s_its",1);
  int max_nda_its   = _my_input_file("Execution/max_nda_its",1);
  
  //create pointer to eigenvector
  Teuchos::RCP<Epetra_Vector> soln = _problem.getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  double keff_err = 1.0;
  double keff_old = keff;
  double keff_older = keff_old;
  
  double keff_err_outer = 1.0;
  double keff_old_outer = keff;
  
  ////--------------------------------------------------
  // first normalize fission source (1/k)*nusigf*phi = 1
  double src_val = _problem.computeFissionSourceStrength(keff, &(*_problem.getSolution()));
  src_val = 1.0/src_val;
  
  _problem.getSolution()->Scale(src_val);
  _problem.copyOldSolution((*_problem.getSolution()));

  for(unsigned int power_its=0;power_its<max_power_its;++power_its)
  {
    std::cout << "power_its: " << power_its<< " " << keff_err<<std::endl;
    keff_older = keff_old;
    keff_old = keff;
    _problem.setKeff(keff);
    
    ////----------------
    // solve LO problem.
    _solver->step();

    ////---------------
    // compute new keff
    const NOX::Epetra::Group& intermediateGroup = dynamic_cast<const NOX::Epetra::Group&>(_solver->getSolutionGroup());
    const Epetra_Vector& intermediateSolution = (dynamic_cast<const NOX::Epetra::Vector&>(intermediateGroup.getX())).getEpetraVector();
    
    keff = _problem.computeKeff(&intermediateSolution);

    ////------------
    // normlize soln
    (*soln) = intermediateSolution;
    double src_val = _problem.computeFissionSourceStrength(keff, &(*soln));
    std::cout<<"in powerMethod srcval "<<src_val<<std::endl;
    //soln->Print(std::cout);

    src_val = 1.0/src_val;
    soln->Scale(src_val);
    
    _problem.copyOldSolution(*soln);

    _solver->reset(noxSoln);
    
    ////----------------
    // check convergence
    keff_err = fabs(keff-keff_old)/keff;

    std::cout.precision(10);
    std::cout << "keff: " << keff << "  error= "<< keff_err<<std::endl;

    if( power_its > 5 )
    {
      double dominance_ratio = (keff-keff_old)/(keff_old-keff_older);
      std::cout << "dominance ratio: " << dominance_ratio <<std::endl;
    }
    
    if( keff_err < keff_tol )
      break;

  }//power_its
  
  return true;
  
}



