
# include "PInterForce_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PInterForce_HertzContact.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PInterForce* PInterForce::PInterForceFactory(const Params& p)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string fij_type = p.input_params("BDApp/inter_particle_forces/type","hertzContact");

   if (fij_type == "hertzContact") return new HertzContact(p);

}
