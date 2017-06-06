
# ifndef FFTSIMPLEDOUBLEWELL_H
# define FFTSIMPLEDOUBLEWELL_H

# include "FFTChemPot_Interface.hpp"

/*
   This is an implementation of the 'FFTChemPot' interface class.
   A double-well function for chemical energy: w*(c^2*(1-c)^2).
*/

class FFTSimpleDoubleWell: public FFTChemPot {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   double w;
   double kap;
   FFTArrays const* const* c;
   const vector<double>& eta;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTSimpleDoubleWell(const Params& p, FFTArrays const* const* cin,
                       const vector<double>& etain) :
                       eta(etain)
   {
      c = cin;
      nc = p.nc;
      w = p.input_params("CHFFTApp/free_energy/w",1.0);
      kap = p.input_params("CHFFTApp/free_energy/kap",1.0);
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTSimpleDoubleWell()
   {
   }

   // -----------------------------
   // Function to calculate mu...
   // -----------------------------

   double muFunc(int n, int ndx, int i, int j, int k)
   {
      double cc = c[n]->getConc(ndx);
      return w*(4*cc*cc*cc - 6*cc*cc + 2*cc);
   }

};

# endif  // FFTSIMPLEDOUBLEWELL_H
