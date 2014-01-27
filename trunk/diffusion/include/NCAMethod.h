#ifndef _NCAMETHOD_H
#define _NCAMETHOD_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "NOX.H"
#include "NOX_Epetra.H"

#include "IterativeScheme.h"

// Finite Volume Problem Class
class NCAMethod : public IterativeScheme
{ 

public:

  /*
   * Constractor with inputfile
   */
  NCAMethod(TwoDProblem &problem,Teuchos::RCP<NOX::Solver::Generic> solver, InputOutput &input_output);
  
  /*
   * Destructor
   */
  ~NCAMethod();


  /*
   * initializing the power iteration.
   */
  void init();

  /*
   * setting up the power iteration.
   */
  void setup();

  /*
   * run the power iteration.
   */
  bool run();
  
private:
  
};


#endif
