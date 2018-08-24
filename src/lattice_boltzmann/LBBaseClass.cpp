
# include "LBBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "LBTypes/scsp2D.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

LBBaseClass* LBBaseClass::LBFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

    // -----------------------------------
    // identify the requested object:
    // -----------------------------------

    string lb_type = input_params("LBApp/type","scsp2D");

    // -----------------------------------
    // return the requested object:
    // -----------------------------------

    if (lb_type == "scsp2D") return new scsp2D(p,input_params);

}
