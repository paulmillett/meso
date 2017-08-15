
# include "PFBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PFTypes/CHBasic.hpp"
# include "PFTypes/CHBD.hpp"
# include "PFTypes/CHBDThinFilm.hpp"
# include "PFTypes/TIPS.hpp"
# include "PFTypes/TIPS2.hpp"
# include "PFTypes/TIPS3.hpp"
# include "PFTypes/TIPS4.hpp"
# include "PFTypes/BCPZone.hpp"
# include "PFTypes/CHTernary.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PFBaseClass* PFBaseClass::PFFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

    // -----------------------------------
    // identify the requested object:
    // -----------------------------------

    string pf_type = input_params("PFApp/type","CHBasic");

    // -----------------------------------
    // return the requested object:
    // -----------------------------------

    if (pf_type == "CHBasic") return new CHBasic(p,input_params);
    if (pf_type == "CHBD")    return new CHBD(p,input_params);
    if (pf_type == "CHBDThinFilm")    return new CHBDThinFilm(p,input_params);
    if (pf_type == "TIPS")    return new TIPS(p,input_params);
    if (pf_type == "TIPS2")    return new TIPS2(p,input_params);
    if (pf_type == "TIPS3")    return new TIPS3(p,input_params);
    if (pf_type == "TIPS4")    return new TIPS4(p,input_params);
    if (pf_type == "BCPZone") return new BCPZone(p,input_params);
    if (pf_type == "CHTernary") return new CHTernary(p,input_params);

}
