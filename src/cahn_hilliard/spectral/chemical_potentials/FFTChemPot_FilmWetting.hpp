
# ifndef FFTFILMWETTING_H
# define FFTFILMWETTING_H

# include "FFTChemPot_Interface.hpp"

/*
   This is an implementation of the 'FFTChemPot' interface class.
   A double-well function with bias on the top- and bottom- z-surfaces.
*/

class FFTFilmWetting: public FFTChemPot {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int nz;
   double w;
   double kap;
   double bias_str;
   FFTArrays const* const* c;
   const vector<double>& eta;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTFilmWetting(const Params& p, FFTArrays const* const* cin,
                  const vector<double>& etain) :
                  eta(etain)
   {
      c = cin;
      nc = p.nc;
      nz = p.nz;
      w = p.input_params("CHFFTApp/free_energy/w",1.0);
      kap = p.input_params("CHFFTApp/free_energy/kap",1.0);
      bias_str = p.input_params("CHFFTApp/free_energy/bias_str",1.0);
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTFilmWetting()
   {
   }

   // -----------------------------
   // Function to calculate mu...
   // -----------------------------

   double muFunc(int n, int ndx, int i, int j, int k)
   {
      double cc = c[n]->getConc(ndx);
      double chempot = w*(cc*cc*cc - cc);
      double bias = 0.0;
      if (k <= 4) bias = double(4-k)/4.0;
      if (k >= (nz-5)) bias = double(k - (nz-5))/4.0;
      chempot += bias_str*2.0*bias*w*(cc+1.0);
      return chempot;
   }

};

# endif  // FFTFILMWETTING_H
