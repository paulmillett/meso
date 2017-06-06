
# include "FFTInitCond_Interface.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "FFTInitCond_HomogeneousRandomNoise.hpp"
# include "FFTInitCond_HomogeneousRandomNoise2.hpp"
# include "FFTInitCond_EmulsionInColloidalCrystal.hpp"
# include "FFTInitCond_SmoothDoubleInterface.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

FFTInitCond* FFTInitCond::FFTInitCondFactory(const Params& p,
                                             FFTArrays** c,
                                             vector<double>& eta)
{

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   string ic_type = p.input_params("CHFFTApp/initial_condition/type","homogeneousRandomNoise");

   if (ic_type == "homogeneousRandomNoise") return new FFTHomogeneousRandomNoise(p,c,eta);
   if (ic_type == "homogeneousRandomNoise2") return new FFTHomogeneousRandomNoise2(p,c,eta);
   if (ic_type == "emulsionInColloidalCrystal") return new FFTEmulsionInColloidalCrystal(p,c,eta);
   if (ic_type == "smoothDoubleInterface") return new FFTSmoothDoubleInterface(p,c,eta);

}
