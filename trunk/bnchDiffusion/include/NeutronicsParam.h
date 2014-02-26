#ifndef NEUTRONICSPARAM_H
#define NEUTRONICSPARAM_H

#include <string>
////
// problem specific includes

////
// generic includes

namespace NeutronicsParam
{
  
  ////
  // global variables
  /**
   * number of groups.
   */
  extern int n_prompt_groups;
  extern int n_precur_groups;

  extern double keff;

  extern std::string execution_type;
  
  /**
   * flag for material reinit
   */
  extern bool need_mat_reinit;
  extern bool eval_ang_flux;
  

}

#endif // NEUTRONICSPARAM_H
