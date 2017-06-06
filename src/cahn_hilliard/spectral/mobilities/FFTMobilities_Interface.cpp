
# include "FFTMobilities_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTMobilities_ZoneAnnealing.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTMobilities* FFTMobilities::FFTMobilitiesFactory(const Params& p, FFTArrays const* const* c)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mob_type = p.input_params("CHFFTApp/mobility/type","isotropic");

   if (mob_type == "zoneAnnealing") return new FFTZoneAnnealing(p,c);

}
