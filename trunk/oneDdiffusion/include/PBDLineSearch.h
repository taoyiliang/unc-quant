#ifndef PBPLINESEARCH_H
#define PBPLINESEARCH_H

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"

#include "NOX_LineSearch_Generic.H"
#include "NOX_Abstract_Group.H"
#include "NOX_LineSearch_FullStep.H"

//forward declaration


class PBDLineSearch : public NOX::LineSearch::Generic
{
public:  
  // constractor that takes the standard line search arguments.
  PBDLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd, 
                Teuchos::ParameterList& params);

  //Destractor
  virtual ~PBDLineSearch();
  
  bool compute(NOX::Abstract::Group& grp,
          double &step,
          const NOX::Abstract::Vector &dir,
          const NOX::Solver::Generic &s);

private:
  NOX::LineSearch::FullStep * _full_step;
  
    
};
#endif
