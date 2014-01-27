#include "PBDLineSearch.h"
#include "NeutronicsParam.h"

#include "NOX.H"
#include "NOX_Epetra.H"
// constractor that takes the standard line search arguments.
PBDLineSearch::PBDLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd, 
                             Teuchos::ParameterList& params)
{
  std::cout << "inside PBDLineSearch:: constractor:" <<std::endl;
  
  //start  with empty constractor.
  _full_step = new NOX::LineSearch::FullStep(gd, 
                                             params);
}


//Destractor
PBDLineSearch:: ~PBDLineSearch()
{
  // empty
}
  
bool PBDLineSearch::compute(NOX::Abstract::Group& grp,
                            double &step,
                            const NOX::Abstract::Vector &dir,
                            const NOX::Solver::Generic &s)
{
//_full_step->compute(grp,step,dir,s);
  
  std::cout << "inside PBDLineSearch::compute "<<step<<std::endl;

  const NOX::Abstract::Group& oldGrp = s.getPreviousSolutionGroup();

  int n_elem = dir.length();
  double max_fraction = 0.75;
  double facmin = 1.0;

  const Epetra_Vector& x = (dynamic_cast<const NOX::Epetra::Vector&>(oldGrp.getX())).getEpetraVector();
  const Epetra_Vector& y = (dynamic_cast<const NOX::Epetra::Vector&>(dir)).getEpetraVector();

  int nelm = x.Map().NumMyElements();
  
  for(unsigned int i = 0;i<nelm;++i)
  {
    double w = x[i]+y[i];

    
//    if( w < (1-max_fraction)*x[i] )
    if( w < 0 )
    {
//      double fac = -max_fraction*x[i]/y[i];
      if( y[i] !=0 && x[i] !=0)
      {
        double fac = -max_fraction*x[i]/y[i];
        if( fac < facmin)
          facmin = fac;
//        std::cout << "i: "<< i << " " << w << " " << x[i]<< " " << y[i]<< " " << facmin<< " " << fac <<std::endl;
      }
      
    }
//    else if( w > (1+max_fraction)*x[i] )
//    {
//      double fac = max_fraction*x[i]/y[i];
//      if( fac < facmin)
    //       facmin = fac;
//    }
  }

  facmin=1;
  
  std::cout << "facmin: " << facmin << std::endl;

  step = facmin;
  grp.computeX(oldGrp,dir,step);
//  NeutronicsParam::eval_ang_flux = true;
  
  return true;
}

    
