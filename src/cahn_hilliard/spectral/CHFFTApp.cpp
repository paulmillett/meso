
# include "CHFFTApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHFFTApp::CHFFTApp(const GetPot& input_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

   p.input_params = input_params;
   p.nc = input_params("CHFFTApp/nc",1);
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

   //	---------------------------------------
   //	Use FFTw to decide 'nx' & 'xOff':
   //	---------------------------------------

   ptrdiff_t locsize, locnx, offx;
   int gxyz = fftw_mpi_local_size_3d(p.NX,p.NY,p.NZ/2+1,MPI_COMM_WORLD,&locnx,&offx);
   p.nx = locnx;
   p.xOff = offx;

   //	---------------------------------------
   //	Create a CHSpecies object:
   //	---------------------------------------

   ch_system = new CHFFTSystem(p);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHFFTApp::~CHFFTApp()
{
   delete ch_system;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void CHFFTApp::initSystem()
{
   ch_system->initializeCahnHilliard();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void CHFFTApp::stepForward(int step)
{

   // ----------------------------------------
   //	Set the time step:
   // ----------------------------------------

   current_step = step;
   ch_system->setTimeStep(current_step);

   // ----------------------------------------
   //	Update Cahn-Hilliard:
   // ----------------------------------------

   ch_system->updateCahnHilliard();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void CHFFTApp::writeOutput(int step)
{
   ch_system->writeOutputFiles(step);
}
