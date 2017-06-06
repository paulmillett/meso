
# include "FFTFilmWetting.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FFTFilmWetting::FFTFilmWetting(const Params& p, FFTArrays** cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   c = cin;
   nc = p.nc;
   nz = p.nz;
   w = p.input_params("CHFFTApp/free_energy/w",1.0);
   kap = p.input_params("CHFFTApp/free_energy/kap",1.0);
   bias_str = p.input_params("CHFFTApp/free_energy/bias_str",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FFTFilmWetting::~FFTFilmWetting()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

void FFTFilmWetting::muFunc(int n, int ndx, int i, int j, int k)
{
   double cc = c[n]->getConc(ndx);
   double chempot = w*(cc*cc*cc - cc);
   double bias = 0.0;
   if (k <= 4) bias = double(4-k)/4.0;
   if (k >= (nz-5)) bias = double(k - (nz-5))/4.0;
   chempot += bias_str*2.0*bias*w*(cc+1.0);
   c[n]->setMu(ndx,chempot);
}
