
# include "FFTChemPot_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTChemPot_SimpleDoubleWell.hpp"
# include "FFTChemPot_FilmWetting.hpp"
# include "FFTChemPot_EmulsionWettingColloids.hpp"
# include "FFTChemPot_TripleWell.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTChemPot* FFTChemPot::FFTChemPotFactory(const Params& p,
                                          FFTArrays const* const* c,
                                          const vector<double>& eta)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string mu_type = p.input_params("CHFFTApp/free_energy/type","simpleDoubleWell");

   if (mu_type == "simpleDoubleWell") return new FFTSimpleDoubleWell(p,c,eta);
   if (mu_type == "filmWetting") return new FFTFilmWetting(p,c,eta);
   if (mu_type == "emulsionWettingColloids") return new FFTEmulsionWettingColloids(p,c,eta);
   if (mu_type == "tripleWell") return new FFTTripleWell(p,c,eta);

}
