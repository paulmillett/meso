
# include "PInitCond_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PInitCond_CubicLattice.hpp"
# include "PInitCond_Random.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PInitCond* PInitCond::PInitCondFactory(const GetPot& p,
                                       vector<double>& r,
                                       vector<double>& v,
                                       vector<double>& rad)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string ic_type = p("BDApp/initial_condition/type","cubicLattice");

   if (ic_type == "cubicLattice") return new CubicLattice(p,r,v,rad);
   if (ic_type == "random") return new Random(p,r,v,rad);

}
