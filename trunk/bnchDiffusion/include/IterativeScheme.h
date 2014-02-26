/*
 * Class Iteative Scheme 
 * This class is a virtual class, which control various eigenvalue iteration, 
 *
 */ 

#ifndef _ITERATIVESCHEME_H
#define _ITERATIVESCHEME_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "NOX.H"
#include "NOX_Epetra.H"

#include "TwoDProblem.h"
#include "InputOutput.h"

#include "getpot.h"
// Finite Volume Problem Class
class IterativeScheme
{ 

public:

  /*
   * Constractor with inputfile
   */
  IterativeScheme(TwoDProblem &problem,Teuchos::RCP<NOX::Solver::Generic> solver, InputOutput &input_output);
  
  /*
   * Destructor
   */
  ~IterativeScheme();


  /*
   * initializing the power iteration.
   */
  virtual void init();

  /*
   * setting up the power iteration.
   */
  virtual void setup();

  /*
   * run the power iteration.
   */
  virtual bool run() =0;

  /*
   * checking the convergence in keff.
   */
  virtual bool checkKeffConvergence();


  /*
   * reinitializing solution vector
   */
  virtual void reinit();

  /*
   * setting NOX output utilities.
   */ 

  
  virtual void setNOXUtils(NOX::Utils &utils);
  

protected:
  TwoDProblem & _problem;

  Teuchos::RCP<NOX::Solver::Generic> _solver;

  InputOutput &_input_output;
  
  NOX::Utils *_utils;

  GetPot &_my_input_file;
  
};


#endif




