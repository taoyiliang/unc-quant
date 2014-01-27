#ifndef MAPS_H
#define MAPS_H

////
// forward declarations.
class Epetra_Map;

class Maps
{
public:
  /**
   * constractor
   */
  Maps(Epetra_Map &standardMap,
       Epetra_Map &overlapMap,
       Epetra_Map &standardDomainMap,
       Epetra_Map &overlapDomainMap,
       std::vector<Epetra_Map *> &standardEqnMap,
       std::vector<Epetra_Map *> &overlapEqnMap);
  /**
   * destractor
   */
  virtual ~Maps()
  {};

  /**
   * this class generally don't do anything but holds all the maps created
   * so that we can access easily
   */
  ////
  // map contains overall dofs. (spans 1~ ncells*neq)
  Epetra_Map &_standardMap; 
  Epetra_Map &_overlapMap;

  ////
  // map contains domain decomposion info using the first pde. (spans 1~ncells)
  Epetra_Map &_standardDomainMap; 
  Epetra_Map &_overlapDomainMap;

  ////
  //vector of map  for cell number to overall dof number.
  std::vector<Epetra_Map *> _standardEqnMap;
  std::vector<Epetra_Map *> _overlapEqnMap;

};

#endif
