
# include "CHApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHApp::CHApp(const GetPot& input_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

   p.input_params = input_params;
   p.nc = input_params("CHApp/nc",1);
   p.NX = input_params("Domain/nx",1);
   p.NY = input_params("Domain/ny",1);
   p.NZ = input_params("Domain/nz",1);
   p.dx = input_params("Domain/dx",1.0);
   p.dy = input_params("Domain/dy",1.0);
   p.dz = input_params("Domain/dz",1.0);
   p.dt = input_params("Time/dt",1.0);
   p.ny = p.NY;
   p.nz = p.NZ;

   //	---------------------------------------
   //	Create cartesian grid of processors:
   //	---------------------------------------

   p.np = MPI::COMM_WORLD.Get_size();   // # of processors
   p.rank = MPI::COMM_WORLD.Get_rank(); // my processor number
   p.nbrL = (p.rank-1) + ((p.rank-1) < 0)*p.np;        // left proc. neighbor
   p.nbrR = (p.rank+1) - ((p.rank+1) > (p.np-1))*p.np; // right proc. neighbor

   //	---------------------------------------
   //	Calculate my sub-grid dimensions:
   //	---------------------------------------

   int nxAverage = int(p.NX/p.np);
	int nxExtra = p.NX%p.np;

   // my 'nx' portion of grid...
	if ((p.rank+1) <= nxExtra) p.nx = nxAverage + 1;
	if ((p.rank+1)  > nxExtra) p.nx = nxAverage;

   // my x-offset...
	int offset = 1;
	for (int i=0; i<p.rank; i++) {
		if (i+1 <= nxExtra) offset += nxAverage + 1;
		if (i+1  > nxExtra) offset += nxAverage;
	}
	p.xOff = offset;

   //	---------------------------------------
   //	Create a CHSpecies object:
   //	---------------------------------------

   ch_system = new CHSystem(p);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHApp::~CHApp()
{
   delete ch_system;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void CHApp::initSystem()
{
   ch_system->initializeCahnHilliard();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void CHApp::stepForward(int step)
{

   // ----------------------------------------
   //	Set the time step:
   // ----------------------------------------

   current_step = step;
   ch_system->setTimeStep(current_step);

   // ----------------------------------------
   //	Update CH system:
   // ----------------------------------------

   ch_system->updateCahnHilliard();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void CHApp::writeOutput(int step)
{
   ch_system->writeOutputFiles(step);
}
