
# ifndef FFTTRIPLEWELL_H
# define FFTTRIPLEWELL_H

# include "FFTChemPot_Interface.hpp"

/*
   This is an implementation of the 'FFTChemPot' interface class.
   A triple-well function for A-B-C mixtures, where A and C fluids do not touch.
*/

class FFTTripleWell: public FFTChemPot {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int nz;
   double w;
   double wab;
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

   FFTTripleWell(const Params& p, FFTArrays const* const* cin,
                 const vector<double>& etain) :
                 eta(etain)
   {
      c = cin;
      nc = p.nc;
      nz = p.nz;
      w = p.input_params("CHFFTApp/free_energy/w",1.0);
      wab = p.input_params("CHFFTApp/free_energy/wab",1.0);
      kap = p.input_params("CHFFTApp/free_energy/kap",1.0);
      bias_str = p.input_params("CHFFTApp/free_energy/bias_str",0.0);
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTTripleWell()
   {
   }

   // -----------------------------
   // Function to calculate mu...
   // -----------------------------

   double muFunc(int n, int ndx, int i, int j, int k)
   {
      int m = 1 - n;
      double cn = c[n]->getConc(ndx);
      double cm = c[m]->getConc(ndx);
      double chempot = w*(4*cn*cn*cn - 6*cn*cn + 2*cn + 2*wab*cn*cm*cm);
      double bias = 0.0;
      if (k <= 4) bias = double(4-k)/4.0;
      if (k >= (nz-5)) bias = double(k - (nz-5))/4.0;
      chempot += bias_str*2.0*bias*w*(cn);
      return chempot;
   }

};

# endif  // FFTTRIPLEWELL_H
