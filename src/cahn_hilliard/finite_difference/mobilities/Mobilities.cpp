
# include "Mobilities.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "ZoneAnnealing.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

Mobilities* Mobilities::MobilitiesFactory(const Params& p,
                                          const vector<FDVector>& c)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mob_type = p.input_params("CHApp/mobility/type","none");

   if (mob_type == "zoneAnnealing") return new ZoneAnnealing(p,c);

}
