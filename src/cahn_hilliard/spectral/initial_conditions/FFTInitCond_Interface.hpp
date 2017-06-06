
# ifndef FFTINITCOND_H
# define FFTINITCOND_H

# include <vector>
# include "../FFTArrays.hpp"
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for source terms in the FFTCH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class FFTInitCond {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of FFTSourceTerms
   // sub-classes:
   // -------------------------------------------------------------------

   static FFTInitCond* FFTInitCondFactory(const Params&,
                                          FFTArrays**,
                                          vector<double>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual void icFunc() = 0;

};

# endif  // FFTINITCOND_H
