

# ifndef FFTDIBLOCKCOPOLYMER_H
# define FFTDIBLOCKCOPOLYMER_H

# include "FFTSourceTerms_Interface.hpp"

/*
   This is an implementation of the 'FFTSourceTerms' interface class.
   This class defines a function that adds a block-copolymer source term to
   the C-H equation to penalize domain growth.
*/

class FFTDiblockCopolymer: public FFTSourceTerms {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   double dt;
   double A0;
   double c0_init_mean;
   FFTArrays const* const* c;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTDiblockCopolymer(const Params& p, FFTArrays const* const* cin)
   {
      c = cin;
      nc = p.nc;
      dt = p.dt;
      A0 = p.input_params("CHFFTApp/source/A0",1.0);
      c0_init_mean = p.input_params("CHFFTApp/initial_condition/c0_init_mean",0.5);
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTDiblockCopolymer()
   {
   }

   // -----------------------------
   // Function to calculate s.t.
   // -----------------------------

   double srcFunc(int n,int i)
   {
      double cc = c[n]->getConc(i);
      return -A0*(cc - c0_init_mean);
   }

};

# endif  // FFTDIBLOCKCOPOLYMER_H
