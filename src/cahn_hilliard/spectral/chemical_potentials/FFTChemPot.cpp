
# include "FFTChemPot.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTSimpleDoubleWell.hpp"
# include "FFTFilmWetting.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTChemPot* FFTChemPot::FFTChemPotFactory(const Params& p, FFTArrays* c[])
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mu_type = p.input_params("CHFFTApp/free_energy/type","simpleDoubleWell");

   if (mu_type == "simpleDoubleWell") return new FFTSimpleDoubleWell(p,c);
   if (mu_type == "filmWetting") return new FFTFilmWetting(p,c);

}
