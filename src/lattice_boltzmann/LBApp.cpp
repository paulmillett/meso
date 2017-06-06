
# include "LBApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

LBApp::LBApp(const GetPot& InputParams)
{

   //	---------------------------------------
   //	create a CHSpecies object:
   //	---------------------------------------

   cout << "Hello from LBApp..." << endl;

   //lb_system = new LBSystem(InputParams);
   mcmp_object = new mcmp(InputParams);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBApp::~LBApp()
{
   //delete lb_system;
   delete mcmp_object;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void LBApp::initSystem()
{
   current_step = 0;
   int iskip = 1;
   int jskip = 1;
   int kskip = 1;
   //lb_system->parseInitialCondition();
   //lb_system->writeVTKFile("rho",current_step,iskip,jskip,kskip);
   mcmp_object->initializeMCMP();
   mcmp_object->writeVTKFile("rhoA",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void LBApp::stepForward(int step)
{

   current_step = step;
   mcmp_object->setTimeStep(current_step);
   mcmp_object->updateMCMP();

   // // ----------------------------------------
   // //	Set the time step:
   // // ----------------------------------------
   //
   // current_step = step;
   // lb_system->setTimeStep(current_step);
   //
   // // ----------------------------------------
   // //	Calculate macro's:
   // // ----------------------------------------
   //
   // lb_system->calculateMacros();
   //
   // // ----------------------------------------
   // //	Calculate fluid forces:
   // // ----------------------------------------
   //
   // lb_system->parseFluidForces();
   //
   // // ----------------------------------------
   // //	Equilibrium distribution:
   // // ----------------------------------------
   //
   // lb_system->equilibriumDistribution();
   //
   // // ----------------------------------------
   // //	Collision:
   // // ----------------------------------------
   //
   // lb_system->collisionStep();
   //
   // // ----------------------------------------
   // //	Streaming:
   // // ----------------------------------------
   //
   // lb_system->streamingStep();
   //
   // // ----------------------------------------
   // //	Bounce-back off solid surfaces:
   // // ----------------------------------------
   //
   // lb_system->bounceBack();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void LBApp::writeOutput(int step)
{
   int iskip = 1;
   int jskip = 1;
   int kskip = 1;
   //lb_system->writeVTKFile("rho",step,iskip,jskip,kskip);
   mcmp_object->writeVTKFile("rhoA",current_step,iskip,jskip,kskip);
}
