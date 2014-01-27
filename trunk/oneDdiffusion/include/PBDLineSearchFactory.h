#ifndef PBDLINESEARCHFACTORY_H
#define PBDLINESEARCHFACTORY_H

#include "NOX_LineSearch_UserDefinedFactory.H"


class PBDLineSearchFactory : public NOX::LineSearch::UserDefinedFactory
{
public:
  PBDLineSearchFactory();
  
  ~PBDLineSearchFactory();
  
  Teuchos::RCP<NOX::LineSearch::Generic> buildLineSearch(const Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params) const;

private:
  Teuchos::RCP<NOX::LineSearch::Generic> _pbd_line_search;
  
  
    
};
#endif
