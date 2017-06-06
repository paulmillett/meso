
# ifndef FFTCHEMPOT_H
# define FFTCHEMPOT_H

# include <vector>
# include "../FFTArrays.hpp"
# include "../Params.hpp"

// ---------------------------------------------------------------------
// This is the base class for chemical potentials in the CH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class FFTChemPot {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of ChemPot sub-classes:
   // -------------------------------------------------------------------

   static FFTChemPot* FFTChemPotFactory(const Params&,
                                        FFTArrays const* const*,
                                        const vector<double>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double muFunc(int,int,int,int,int) = 0;

};

# endif  // FFTCHEMPOT_H
