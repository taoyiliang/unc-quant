////
//trilino  includes
#include "Epetra_Map.h"

////
//generic includes
#include <iostream>
#include <vector>

#include "Maps.h"

Maps::Maps(Epetra_Map &standardMap,
           Epetra_Map &overlapMap,
           Epetra_Map &standardDomainMap,
           Epetra_Map &overlapDomainMap,
           std::vector<Epetra_Map *> &standardEqnMap,
           std::vector<Epetra_Map *> &overlapEqnMap)
  :_standardMap(standardMap), 
   _overlapMap(overlapMap),
   _standardDomainMap(standardDomainMap),
   _overlapDomainMap(overlapDomainMap)
{
  int nmap = standardEqnMap.size();
  _standardEqnMap.resize(nmap);
  _overlapEqnMap.resize(nmap);
  
  for(unsigned int imap=0;imap<nmap;++imap)
  {
    _standardEqnMap[imap] = standardEqnMap[imap];
    _overlapEqnMap[imap] = overlapEqnMap[imap];
  }
}
