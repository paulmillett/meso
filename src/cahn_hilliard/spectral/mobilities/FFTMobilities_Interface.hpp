
# ifndef FFTMOBILITIES_H
# define FFTMOBILITIES_H

# include "../FFTArrays.hpp"
# include "../Params.hpp"

// ---------------------------------------------------------------------
// This is the base class for chemical potentials in the CH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class FFTMobilities {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of ChemPot sub-classes:
   // -------------------------------------------------------------------

   static FFTMobilities* FFTMobilitiesFactory(const Params&, FFTArrays const* const*);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double mobFunc(int,int,int,int,int) = 0;

};

# endif  // FFTMOBILITIES_H
