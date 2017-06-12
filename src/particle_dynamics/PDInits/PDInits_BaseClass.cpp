
# include "PDInits_BaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PDInits_CubicLattice.hpp"
# include "PDInits_Random.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PDInits_BaseClass* PDInits_BaseClass::PDInitFactory(const GetPot& p,
                                                    vector<double>& r,
                                                    vector<double>& v,
                                                    vector<double>& rad)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string ic_type = p("PDApp/initial_condition/type","cubicLattice");

   if (ic_type == "cubicLattice") return new CubicLattice(p,r,v,rad);
   if (ic_type == "random") return new Random(p,r,v,rad);

}
