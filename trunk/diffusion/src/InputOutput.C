////
// generic includes
#include <iostream>

////
// problem specific includes
#include "FVHelper.h"
#include "InputOutput.h"
#include "NeutronicsParam.h"

////
// trilinos includes
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

////
// local
#include "TwoDMesh.h"

InputOutput::InputOutput(GetPot &input_file)
  :_my_input_file(input_file)
{}


void InputOutput::write2File(std::string outfile, Epetra_Vector& soln, TwoDMesh &mesh,double k)
{
  std::ofstream my_output_file;
  my_output_file.precision(10);
  my_output_file<<std::scientific;
  my_output_file.open (outfile.c_str());

  my_output_file << "DIM,"<<mesh.nX<<","<<mesh.nY<<std::endl;
  my_output_file << "k,"<<k<<std::endl;

  my_output_file << "m1siga2,"<<mesh.materials[1].siga[1]<<std::endl;

  for( int icell=0;icell<soln.MyLength()/2;++icell)
  {
    std::vector<int> pos=mesh.ij(icell);
    int gid1=mesh.gid(pos[0],pos[1],0);
    int gid2=mesh.gid(pos[0],pos[1],1);
    my_output_file << mesh.cells[icell].pos[0]<<","<<mesh.cells[icell].pos[1] <<" | " <<  soln[gid1]<<"  , "<<soln[gid2]<<std::endl;
  }

  return;
  
}

