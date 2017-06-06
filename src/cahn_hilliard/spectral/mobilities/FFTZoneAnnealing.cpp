
# include "FFTZoneAnnealing.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FFTZoneAnnealing::FFTZoneAnnealing(const Params& p, FFTArrays** cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   c = cin;
   nc = p.nc;
   vzone = p.input_params("CHFFTApp/mobility/vzone",1.0);
   wzone = p.input_params("CHFFTApp/mobility/wzone",20.0);
   
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FFTZoneAnnealing::~FFTZoneAnnealing()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

void FFTZoneAnnealing::mobFunc(int n, int ndx, int i, int j, int k)
{
   // need to define later
}
