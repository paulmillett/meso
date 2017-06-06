
# ifndef CHEMPOT_H
# define CHEMPOT_H

# include <vector>
# include "../FDVector.hpp"
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for chemical potentials in the CH App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class ChemPot {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of ChemPot sub-classes:
   // -------------------------------------------------------------------

   static ChemPot* ChemPotFactory(const Params&,const vector<FDVector>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double muFunc(int,int,int,int,int) = 0;

};

# endif  // CHEMPOT_H
