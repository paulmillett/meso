
# include "FFTDiblockCopolymer.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FFTDiblockCopolymer::FFTDiblockCopolymer(const Params& p, FFTArrays** cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   c = cin;
   nc = p.nc;
   dt = p.dt;
   A0 = p.input_params("CHFFTApp/source/A0",1.0);
   c0_init_mean = p.input_params("CHFFTApp/initial_condition/c0_init_mean",0.5);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FFTDiblockCopolymer::~FFTDiblockCopolymer()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

void FFTDiblockCopolymer::srcFunc(int n,int i)
{
   double cc = c[n]->getConc(i);
   double fterm = -A0*(cc - c0_init_mean);
   c[n]->addtoConc(i,dt*fterm);
}
