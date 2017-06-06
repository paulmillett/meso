
# ifndef FFTSOURCETERMS_H
# define FFTSOURCETERMS_H

# include <vector>
# include "../FFTArrays.hpp"
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for source terms in the FFTCH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class FFTSourceTerms {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of FFTSourceTerms
   // sub-classes:
   // -------------------------------------------------------------------

   static FFTSourceTerms* FFTSourceTermsFactory(const Params&, FFTArrays*[]);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual void srcFunc(int,int) = 0;

};

# endif  // FFTSOURCETERMS_H
