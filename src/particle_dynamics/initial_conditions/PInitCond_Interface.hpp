
# ifndef PINITCOND_H
# define PINITCOND_H

# include <vector>
# include "../../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for initial conditions in the BD App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class PInitCond {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of PInitCond
   // sub-classes:
   // -------------------------------------------------------------------

   static PInitCond* PInitCondFactory(const GetPot&,
                                      vector<double>&,
                                      vector<double>&,
                                      vector<double>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual void icFunc() = 0;

};

# endif  // PINITCOND_H
