
# include "LBBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "LBTypes/mcmp2D.hpp"
# include "LBTypes/mcmp3D.hpp"

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

	string lb_type = input_params("LBApp/type","mcmp2D");

	// -----------------------------------
	// return the requested object:
	// -----------------------------------

	if (lb_type == "mcmp2D") return new mcmp2D(p,input_params);
	if (lb_type == "mcmp3D") return new mcmp3D(p,input_params);
	
	// if input file doesn't have a correct type return a nullptr
	return NULL;

}
