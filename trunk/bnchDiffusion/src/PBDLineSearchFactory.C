////
////nox include
#include "NOX_LineSearch_Generic.H"
#include "NOX_GlobalData.H"

#include "PBDLineSearchFactory.h"
#include "PBDLineSearch.h"



PBDLineSearchFactory::PBDLineSearchFactory()
{
  std::cout << "PBDLineSearchFactory::constractor:"<<std::endl;
  
  // empty constractor
}

  
PBDLineSearchFactory::~PBDLineSearchFactory()
{
  // empty
}

  
Teuchos::RCP<NOX::LineSearch::Generic>
PBDLineSearchFactory::buildLineSearch(const Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params) const
{
  std::cout << "build linesearch"<<std::endl;
  
  Teuchos::RCP<NOX::LineSearch::Generic> _pbd_line_search = Teuchos::rcp(new PBDLineSearch(gd, params));
  return _pbd_line_search;
  
}

