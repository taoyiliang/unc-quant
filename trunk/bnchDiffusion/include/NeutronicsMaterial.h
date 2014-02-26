c#ifndef _NEUTRONICSMATERIAL_H
#define _NEUTRONICSMATERIAL_H

#include "Teuchos_RCP.hpp"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_IntVector;
class Epetra_Vector;
class Epetra_Import;
class Epetra_Operator;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;


// Finite Volume Problem Class
class NeutronicsMaterial
{ 

public:

  NeutronicsMaterial(std::vector<Epetra_Map* > &map, Epetra_IntVector &mat_id);

  void getXS(std::vector<Epetra_Vector *> &u_vec);
  
private:

};

  

#endif//_NEUTRONICSMATERIAL_H
