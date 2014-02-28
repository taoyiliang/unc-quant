#include <time.h>
#include <sstream>
#include <iostream>
// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_LineSearch_UserDefinedFactory.H"
#include "NOX_Utils.H"
// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// User's application specific files 
#include "getpot.h"             // parsing input file
#include "Problem_Interface.h"  // Interface file to NOX
#include "TwoDProblem.h"              
#include "PhysicsBasedPrecon.h"
#include "MLPrecond.h"
#include "NoPreconditioner.h"
#include "FVHelper.h"
#include "InputOutput.h"
#include "TwoDMesh.h"
#include "Materials.h"
#include "Cells.h"
#include "MatProps.h"
#include "PBDLineSearchFactory.h"

#include "PowerMethod.h"
#include "NCAMethod.h"
#include "NeutronicsParam.h"

// Required for reading and writing parameter lists from xml format
// Configure Trilinos with --enable-teuchos-extended
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

//using namespace std;

int main(int argc, char *argv[])
{
  std::cout.precision(10);
  
  int ierr = 0;

  clock_t start_time,stop_time;
  start_time=clock();
  
  
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  ////---------------------------
  // call getpot based inputfile.
  GetPot command_line(argc,argv);
  
  std::string input_filename = "";
  if ( command_line.search("-i") )
    input_filename = command_line.next(input_filename.c_str());
  else
  {
    std::cout<<"Must specify an input file using -i"<<std::endl;
    throw;
  }

  ////-----------------
  // reading input file
  GetPot input_file(input_filename.c_str() );

  ////------------------
  // I/O class.
  InputOutput input_output(input_file);

  ////----------
  // create mesh
  TwoDMesh mesh;

  // create material classes 
  int nMats = input_file.vector_variable_size("Material/names");
  std::string matName;
  int matLabel;
  std::string prefix;
  std::vector<double> D,siga,sigc,nsigf,sigtr,chi;
  D.resize(2);
  siga.resize(2);
  sigc.resize(2);
  nsigf.resize(2);
  sigtr.resize(2);
  chi.resize(2);
  for (int i=0;i<nMats;++i)
  {
    input_file.set_prefix("");
    matName=input_file("Material/names"," ",i);
    matLabel=input_file("Material/labels",-1,i);
    for (int g=0;g<2;++g)
    {
      int gp = g+1;
      std::ostringstream str_g;
      str_g << (gp);
      prefix=std::string("Material/")+matName+std::string("/g")+
                            str_g.str()+std::string("/");
      //cout<<"prefix: "<<prefix<<endl;
      input_file.set_prefix(prefix);
      D[g]     = input_file("D"    ,1.0);
      nsigf[g] = input_file("nsigf",1.0);
      sigc[g]  = input_file("sigc" ,1.0);
      sigtr[g] = input_file("sigtr",1.0);
      chi[g]   = input_file("chi"  ,1.0);
      siga[g]  = nsigf[g]/2.43+sigc[g];
    }
    mesh.add_material(matLabel,matName,D,siga,nsigf,sigtr,chi);
    //cout<<"  D     "<<mesh.materials[matLabel].D[0]<<" "<<D[1]<<endl;
    //cout<<"  siga  "<<mesh.materials[matLabel].siga[0]<<" "<<siga[1]<<endl;
    //cout<<"  nsigf "<<mesh.materials[matLabel].nsigf[0]<<" "<<nsigf[1]<<endl;
    //cout<<"  sigtr "<<mesh.materials[matLabel].sigtr[0]<<" "<<sigtr[1]<<endl;
    //cout<<"  chi   "<<mesh.materials[matLabel].chi[0]<<" "<<chi[1]<<endl;
  }

  //read in regions
  int ictr=0;
  std::vector<int> mat_id;
  input_file.set_prefix("");
  int NX=input_file("Mesh/nXregs",1);
  int NY=input_file("Mesh/nYregs",1);
  //mat_id.resize(NY*NX);
  for (int j=0;j<NY;++j)
  {
    //mat_id[j].resize(NX);
    for (int i=0;i<NX;++i)
    {
      mat_id.push_back(input_file("Mesh/matmap",1,ictr));
      //mat_id[j][i] = input_file("Mesh/matmap",-1,ictr);
      ictr ++;
    }//i
  }//j

  int nx_per_reg = input_file("Mesh/nx_per_reg",1);
  int ny_per_reg = input_file("Mesh/ny_per_reg",1);
  double dx = input_file("Mesh/len_Xside",1.0)/double(nx_per_reg);
  double dy = input_file("Mesh/len_Yside",1.0)/double(ny_per_reg);
  std::cout<<"dx,dy "<<dx<<" "<<dy<<std::endl;
  int nXregs = input_file("Mesh/nXregs",1);
  int nYregs = input_file("Mesh/nYregs",1);

  //double dx;
  //double dy;
  int    ncells_per_Xreg;
  int    ncells_per_Yreg;

  //dx = dx_;
  //dy = dy_;
  ncells_per_Xreg = nx_per_reg;
  ncells_per_Yreg = ny_per_reg;
  
  mesh.init(nXregs, nYregs, dx, dy,
            ncells_per_Xreg,ncells_per_Yreg,
            mat_id);


  ////-------------------
  // create main problem.
  TwoDProblem problem(mesh, Comm);

  ////-------------------
  // set Boundary condition
  int left, right, bottom, top;
  left = input_file("BCs/left",0);
  right =input_file("BCs/right",0);
  top =input_file("BCs/top",0);
  bottom =input_file("BCs/bottom",0);
  
  problem.setBCs(left,right,bottom,top);

  ////--------------------
  // preconditioning object
  PhysicsBasedPrecon  *PBP;
  PBP = new MLPreconditioner("ml_precon",problem);

  ////------------------------------
  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = problem.getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  double &keff = NeutronicsParam::keff;
  std::string &execution_type = NeutronicsParam::execution_type;
  execution_type = input_file("Execution/execution_type","eigenvalue");
  std::cout << "execution_type: " << execution_type<<std::endl;


  ////------------------------------------------------------------------
  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Begin setting Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 10);
  printParams.set("Output Processor", 0);

  printParams.set("Output Information",
                  NOX::Utils::OuterIterationStatusTest + 
                  NOX::Utils::OuterIteration + 
                  NOX::Utils::Warning+
                  NOX::Utils::LinearSolverDetails);
  

  // Create printing utilities
  NOX::Utils utils(printParams);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");

  ////--------------------
  // physics based damping 
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> ls_factory = Teuchos::rcp(new PBDLineSearchFactory,false);
  
  searchParams.set("Method", "User Defined");
  searchParams.set("User Defined Line Search Factory", ls_factory);

  ////--------------------
  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  ////-----------------------------------------------
  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver",           "GMRES");
  lsParams.set("Preconditioner",          "User Defined");
  lsParams.set("Max Iterations",          150);
  lsParams.set("Size of Krylov Subspace", 150);
  lsParams.set("Tolerance",               1e-2); 
  lsParams.set("Max Age Of Prec",         1000); 
  lsParams.set("Output Frequency",        1 );
  
  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NOX_Epetra_Interface
  Teuchos::RCP<Problem_Interface> interface = Teuchos::rcp(new Problem_Interface(problem));

  ////-----------------------------------------------------
  // set pointer to Physics-based preconditioning operator.
  Teuchos::RCP<Epetra_Operator> M = Teuchos::rcp( PBP, false);

  ////-------------------------------------
  // Create the linear system (matrix-free)
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

  
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<Epetra_CrsMatrix> J = problem.GetMatrix() ;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iPrec, M, noxSoln));
//  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac,J, iPrec, M, noxSoln));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln, linSys)); 

  std::cout<<grp<<std::endl;
  ////
  // setting convergence criteria, should be changed to make input params.
  Teuchos::RCP<NOX::StatusTest::Combo> converged;
  {
    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF>    absresid = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    //std::cout<<"got here"<<std::endl;
    // TODO this step is the problem
    Teuchos::RCP<NOX::StatusTest::NormF>    relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-2));
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-9));
    Teuchos::RCP<NOX::StatusTest::NormWRMS>     wrms = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-3, 1.0e-9));
    
    converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    converged->addStatusTest(absresid);
    converged->addStatusTest(relresid);
    converged->addStatusTest(wrms);
    converged->addStatusTest(update);
  }
  
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(100));
  Teuchos::RCP<NOX::StatusTest::FiniteValue>    fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo>       combo = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  Teuchos::RCP<Teuchos::ParameterList> finalParamsPtr = nlParamsPtr;

  // Create the method
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp, combo, finalParamsPtr);


  ////----------------------------------------------
  // Here is main loop for eigenvalue calculation ..

  std::cout << "executing..." <<execution_type<<std::endl;
  if( execution_type == "NCA" )
  {
    NCAMethod nca_method(problem,solver,input_output);
    nca_method.run();
  }
  else
  {
    PowerMethod power_method(problem,solver,input_output);
    power_method.run();
  }

  std::cout<<"Final keff: 1 + "<<keff-1.0<<std::endl;
  ////----------------------
  // write solution to file.
  std::string outfile = input_file("Output/file","test.out");  
  input_output.write2File(outfile, (*soln), mesh,keff);
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif


  stop_time = clock();
  double time = (double)(stop_time-start_time)/CLOCKS_PER_SEC;
  std::cout << " runtime: " << time << std::endl;
  
/* end main
 */
    return ierr ;

}




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


  
