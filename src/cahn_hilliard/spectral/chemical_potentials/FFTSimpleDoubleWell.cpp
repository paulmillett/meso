
# include "FFTSimpleDoubleWell.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FFTSimpleDoubleWell::FFTSimpleDoubleWell(const Params& p, FFTArrays** cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   c = cin;
   nc = p.nc;
   w = p.input_params("CHFFTApp/free_energy/w",1.0);
   kap = p.input_params("CHFFTApp/free_energy/kap",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FFTSimpleDoubleWell::~FFTSimpleDoubleWell()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

void FFTSimpleDoubleWell::muFunc(int n, int ndx, int i, int j, int k)
{
   double cc = c[n]->getConc(ndx);
   double chempot = w*(4*cc*cc*cc - 6*cc*cc + 2*cc);
   c[n]->setMu(ndx,chempot);
}
