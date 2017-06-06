
# include "FFTMobilities.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTZoneAnnealing.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTMobilities* FFTMobilities::FFTMobilitiesFactory(const Params& p, FFTArrays* c[])
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mob_type = p.input_params("CHFFTApp/mobility/type","isotropic");

   if (mob_type == "zoneAnnealing") return new FFTZoneAnnealing(p,c);

}
