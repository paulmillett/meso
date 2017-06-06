
# include "BDApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

BDApp::BDApp(const GetPot& input_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

   p.input_params = input_params;
   p.N  = input_params("BDApp/N",1);
   p.Lx = input_params("BDApp/Lx",5.0);
   p.Ly = input_params("BDApp/Ly",5.0);
   p.Lz = input_params("BDApp/Lz",5.0);
   p.dt = input_params("Time/dt",1.0);

   //	---------------------------------------
   //	Create cartesian grid of processors:
   //	---------------------------------------

   p.np = MPI::COMM_WORLD.Get_size();   // # of processors
   p.rank = MPI::COMM_WORLD.Get_rank(); // my processor number

   //	---------------------------------------
   //	Create a CHSpecies object:
   //	---------------------------------------

   bd_particles = new Particles(p);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

BDApp::~BDApp()
{
   delete bd_particles;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void BDApp::initSystem()
{
   bd_particles->initParticles();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void BDApp::stepForward(int step)
{

   // ----------------------------------------
   //	Set the time step:
   // ----------------------------------------

   current_step = step;
   bd_particles->setTimeStep(current_step);

   // ----------------------------------------
   //	Update Particles:
   // ----------------------------------------

   bd_particles->updateParticles();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void BDApp::writeOutput(int step)
{
   bd_particles->writeVTKFile("particles",step);
}
