#include "IterativeScheme.h"
#include "NeutronicsParam.h"


IterativeScheme::IterativeScheme(TwoDProblem &problem, Teuchos::RCP<NOX::Solver::Generic> solver, InputOutput &input_output) :
  _problem(problem),
  _solver(solver),
  _input_output(input_output),
  _my_input_file(input_output.getInputFile())
{
  
}

IterativeScheme::~IterativeScheme()
{
}

void IterativeScheme::init()
{

}

void IterativeScheme::setup()
{

}

bool IterativeScheme::checkKeffConvergence()
{
  return true;
}

void IterativeScheme::reinit()
{
  Teuchos::RCP<Epetra_Vector> soln = _problem.getSolution();
  _problem.copyOldSolution(*soln);
}



void IterativeScheme::setNOXUtils(NOX::Utils &utils)
{
  _utils = &utils;
}
