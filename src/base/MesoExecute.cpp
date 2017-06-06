
# include "MesoExecute.hpp"
# include "../utils/GetPot"
# include <mpi.h>
# include <cstdlib>
using namespace std;

// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

MesoExecute::MesoExecute()
{

   // -----------------------------------
   // create MPI environment:
   // -----------------------------------

   MPI::Init();
   int np = MPI::COMM_WORLD.Get_size();            // # of processors
   int rank = MPI::COMM_WORLD.Get_rank();          // my processor number

   // -----------------------------------
   // make output directories:
   // -----------------------------------

   if (rank == 0) {
	 	std::system("mkdir -p vtkoutput");            // make output directory
	 	std::system("exec rm -rf vtkoutput/*.vtk");   // remove any existing files
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

MesoExecute::~MesoExecute()
{
   MPI::Finalize();
}



// -------------------------------------------------------------------------
// Create the Meso simulation objects, and store them in a vector:
// -------------------------------------------------------------------------

void MesoExecute::createMesoObjects()
{

   // ------------------------------------------------
   // 'GetPot' object containing input parameters:
   // ------------------------------------------------

   GetPot InParams("input.dat");

   // ------------------------------------------------
   // make vector of input section names:
   // ------------------------------------------------

   vector<string> sections = InParams.get_section_names();

   // ------------------------------------------------
   // determine which sections are executable 'apps':
   // ------------------------------------------------

   for (int i=0; i<sections.size(); i++) {

      // ---------------------------------------------
      // get string that stores value of "section/app"
      // ---------------------------------------------

      string currentSec = sections[i] + "app";
      const char* currentSecChar = currentSec.c_str();
      string appFlag = InParams(currentSecChar,"false");

      // ---------------------------------------------
      // if "app = true", make a new object:
      // ---------------------------------------------

      if (appFlag == "true") {
         mesoapps.push_back(MesoBase::MesoObjectFactory(sections[i]));
      }

   }

   // ------------------------------------------------
   // loop over executable objects, initializing each:
   // ------------------------------------------------

   for (int i=0; i<mesoapps.size(); i++) {
      mesoapps[i]->initSystem();
   }

}



// -------------------------------------------------------------------------
// Execute the simulation by marching forward in time:
// -------------------------------------------------------------------------

void MesoExecute::executeMesoSimulation()
{

   // -----------------------------------
   // get the number of time steps:
   // -----------------------------------

   GetPot InParams("input.dat");
   int nstep = InParams("Time/nstep",0);
   int numOutputs = InParams("Output/numOutputs",1);

   // -----------------------------------
   // determine output interval:
   // -----------------------------------

   int outInterval = 0;
   if (numOutputs != 0) outInterval = nstep/numOutputs;
   if (numOutputs == 0) outInterval = nstep+1;

   // -----------------------------------
   // MARCH THROUGH TIME:
   // -----------------------------------

   for (int step=1; step<nstep+1; step++) {

      // --------------------------------
      // call 'StepForward' for each app:
      // --------------------------------

      for (int i=0; i<mesoapps.size(); i++) {
         mesoapps[i]->stepForward(step);
      }

      // --------------------------------
      // write output for each app:
      // --------------------------------

      if (step == 1 || step%outInterval == 0) {
         for (int i=0; i<mesoapps.size(); i++) {
            mesoapps[i]->writeOutput(step);
         }
      }

   }

}
