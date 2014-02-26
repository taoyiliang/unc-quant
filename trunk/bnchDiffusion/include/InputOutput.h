#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

////
// problem specific includes
#include "NeutronicsParam.h"

////
// generic includes
#include <map>
#include <vector>
#include "getpot.h"

////
// forward declaration
class Epetra_Comm;

class Epetra_Vector;
class TwoDMesh;

class InputOutput
{
 public:
  /**
   * Constractor 
   */
  InputOutput(GetPot &input_file);

  /**
   * Destractor
   */
  virtual ~InputOutput()
  {};
  

  void write2File(std::string outfile,
                  Epetra_Vector &soln,
                  TwoDMesh &mesh,
                  double k=0);
  

  GetPot &getInputFile()
  {
    return _my_input_file;
  }

  
protected:

private:
  GetPot & _my_input_file;
};

#endif //INPUTOUTPUT_H
