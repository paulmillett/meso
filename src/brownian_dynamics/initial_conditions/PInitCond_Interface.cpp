
# include "PInitCond_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PInitCond_CubicLattice.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PInitCond* PInitCond::PInitCondFactory(const Params& p,
                                       vector<double>& r,
                                       vector<double>& v,
                                       vector<double>& rad)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string ic_type = p.input_params("BDApp/initial_condition/type","cubicLattice");

   if (ic_type == "cubicLattice") return new CubicLattice(p,r,v,rad);

}
