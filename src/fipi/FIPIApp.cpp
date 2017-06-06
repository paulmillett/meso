
# include "FIPIApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FIPIApp::FIPIApp(const GetPot& input_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

   p.N  = input_params("FIPIApp/N",1);
   p.NX = input_params("Domain/nx",1);
   p.NY = input_params("Domain/ny",1);
   p.NZ = input_params("Domain/nz",1);
   p.dx = input_params("Domain/dx",1.0);
   p.dy = input_params("Domain/dy",1.0);
   p.dz = input_params("Domain/dz",1.0);
   p.dt = input_params("Time/dt",1.0);
   p.dtRatio = input_params("FIPIApp/dtRatio",1);
   p.visc = input_params("FIPIApp/visc",1.0);
   p.co = input_params("FIPIApp/co",0.0);
   p.pmob = input_params("FIPIApp/pmob",1.0);
   p.rcut = input_params("FIPIApp/rcut",4.0);
   p.A = input_params("FIPIApp/A",1.0);
   p.scl = input_params("FIPIApp/scl",1.0);
   p.Khertz = input_params("FIPIApp/Khertz",1.0);
   p.iskip = input_params("Output/iskip",1);
   p.jskip = input_params("Output/jskip",1);
   p.kskip = input_params("Output/kskip",1);
   p.w = input_params("FIPIApp/w",1.0);
   p.kap = input_params("FIPIApp/kap",1.0);
   p.ny = p.NY;
   p.nz = p.NZ;
   p.Lx = p.NX*p.dx;
   p.Ly = p.NY*p.dy;
   p.Lz = p.NZ*p.dz;

   //	---------------------------------------
   //	Get some MPI parameters:
   //	---------------------------------------

   p.np = MPI::COMM_WORLD.Get_size();   // # of processors
   p.rank = MPI::COMM_WORLD.Get_rank(); // my processor number

   //	---------------------------------------
   //	Use FFTw to decide 'nx' & 'xOff':
   //	---------------------------------------

   ptrdiff_t locsize, locnx, offx;
   fftw_mpi_init();
   locsize = fftw_mpi_local_size_3d(p.NX,p.NY,p.NZ,MPI_COMM_WORLD,&locnx,&offx);
   p.nx = locnx;
   p.xOff = offx;

   //	---------------------------------------
   //	Create a CHSpecies object:
   //	---------------------------------------

   fipi_system = new FIPISystem(p);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FIPIApp::~FIPIApp()
{
   delete fipi_system;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void FIPIApp::initSystem()
{
   fipi_system->initializeFIPI();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void FIPIApp::stepForward(int step)
{

   // ----------------------------------------
   //	Set the time step:
   // ----------------------------------------

   current_step = step;
   fipi_system->setTimeStep(current_step);

   // ----------------------------------------
   //	Update CH system:
   // ----------------------------------------

   fipi_system->updateFIPI();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void FIPIApp::writeOutput(int step)
{
   fipi_system->writeOutputFiles(step);
}
