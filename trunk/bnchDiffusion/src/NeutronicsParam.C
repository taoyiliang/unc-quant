#include "NeutronicsParam.h"

namespace NeutronicsParam
{
  int n_prompt_groups = 0;
  int n_precur_groups = 0;

  double keff = 1.0;
  
  std::string execution_type = "eigenvalue";

  bool need_mat_reinit = false;

  bool eval_ang_flux = true;
    
}
