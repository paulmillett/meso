
# ifndef FFTZONEANNEALING_H
# define FFTZONEANNEALING_H

# include "FFTMobilities_Interface.hpp"

/*
   This is an implementation of the 'FFTMobilities' interface class.
   The mobility term is defined as a moving 'tanh' function.
*/

class FFTZoneAnnealing: public FFTMobilities {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int xOff;
   double dx;
   double vzone;
   double wzone;
   FFTArrays const* const* c;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTZoneAnnealing(const Params& p, FFTArrays const* const* cin)
   {
      c = cin;
      nc = p.nc;
      xOff = p.xOff;
      dx = p.dx;
      vzone = p.input_params("CHFFTApp/mobility/vzone",1.0);
      wzone = p.input_params("CHFFTApp/mobility/wzone",20.0);
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTZoneAnnealing()
   {
   }

   // -----------------------------
   // Function to calculate mob...
   // -----------------------------

   double mobFunc(int n, int ndx, int i, int j, int k)
   {
      double xi = dx*double(i+xOff);
      double xf = 100.0;
      //return 0.5*(1.0 - tanh(6.0*(xi-xf)/wzone));
      return 1.0;
   }

};

# endif  // FFTZONEANNEALING_H
