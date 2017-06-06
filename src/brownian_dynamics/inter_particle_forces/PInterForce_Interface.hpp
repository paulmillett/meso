
# ifndef PINTERFORCE_H
# define PINTERFORCE_H

# include <vector>
# include "../Params.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for inter-particle forces in the BD App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class PInterForce {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of PInterForce
   // sub-classes:
   // -------------------------------------------------------------------

   static PInterForce* PInterForceFactory(const Params&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual double fijFunc(double,double) = 0;

};

# endif  // PINTERFORCE_H
