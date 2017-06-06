

# ifndef HERTZCONTACT_H
# define HERTZCONTACT_H

# include "PInterForce_Interface.hpp"

/*
   This is an implementation of the 'PInterForce' interface class.
   It defines the inter-particle forces based on the Hertz contact model.
*/

class HertzContact : public PInterForce {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   double K;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   HertzContact(const Params& p)
   {
      K = p.input_params("BDApp/inter_particle_forces/K",0.2);
   }


   // -----------------------------
   // Destructor...
   // -----------------------------

   ~HertzContact()
   {
   }


   // -----------------------------
   // Function to calculate fij:
   // -----------------------------

   double fijFunc(double rij, double s2s)
   {
      if (s2s < 0.) {
         return 2.5*K*pow(-s2s,1.5);
      }
      else {
         return 0.0;
      }
   }

};

# endif  // HERTZCONTACT_H
