
# include "SourceTerms.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "None.hpp"
# include "DiblockCopolymer.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

SourceTerms* SourceTerms::SourceTermsFactory(const Params& p,
                                             const vector<FDVector>& c)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string src_type = p.input_params("CHApp/source/type","none");

   if (src_type == "none") return new None();
   if (src_type == "diblockCopolymer") return new DiblockCopolymer(p,c);

}
