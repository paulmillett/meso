
# ifndef SOURCETERMS_H
# define SOURCETERMS_H

# include <vector>
# include "../FDVector.hpp"
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for chemical potentials in the CH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class SourceTerms {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of ChemPot sub-classes:
   // -------------------------------------------------------------------

   static SourceTerms* SourceTermsFactory(const Params&,const vector<FDVector>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double srcFunc(int,int,int,int,int) = 0;

};

# endif  // SOURCETERMS_H
