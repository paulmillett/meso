
# ifndef MOBILITIES_H
# define MOBILITIES_H

# include <vector>
# include "../FDVector.hpp"
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for chemical potentials in the CH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class Mobilities {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of Mobility sub-classes:
   // -------------------------------------------------------------------

   static Mobilities* MobilitiesFactory(const Params&,const vector<FDVector>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double mobFunc(int,int,int,int,int) = 0;

};

# endif  // MOBILITIES_H
