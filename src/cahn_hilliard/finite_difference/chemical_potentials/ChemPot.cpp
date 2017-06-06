
# include "ChemPot.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "SimpleDoubleWell.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

ChemPot* ChemPot::ChemPotFactory(const Params& p,
                                 const vector<FDVector>& c)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mu_type = p.input_params("CHApp/free_energy/type","simpleDoubleWell");

   if (mu_type == "simpleDoubleWell") return new SimpleDoubleWell(p,c);

}
