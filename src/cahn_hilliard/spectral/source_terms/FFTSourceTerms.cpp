
# include "FFTSourceTerms.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTDiblockCopolymer.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTSourceTerms* FFTSourceTerms::FFTSourceTermsFactory(const Params& p,
                                                      FFTArrays* c[])
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string src_type = p.input_params("CHFFTApp/source/type","none");

   if (src_type == "diblockCopolymer") return new FFTDiblockCopolymer(p,c);

}
