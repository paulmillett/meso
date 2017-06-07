
# include "PDBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PDTypes/Hertz.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PDBaseClass* PDBaseClass::PDFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

   // -----------------------------------
   // identify the requested object:
   // -----------------------------------

   string pd_type = input_params("PDApp/type","Hertz");

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   if (pd_type == "Hertz") return new Hertz(p,input_params);

}



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

PDBaseClass::PDBaseClass(const CommonParams& p, const GetPot& input_params)
{
  cout << "Hello from PDBaseClass Ctor" << endl;
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

PDBaseClass::~PDBaseClass()
{

}



// -------------------------------------------------------------------------
// Updater:
// -------------------------------------------------------------------------

void PDBaseClass::updateParticles()
{

}



// -------------------------------------------------------------------------
// Outputer:
// -------------------------------------------------------------------------

void PDBaseClass::outputParticles()
{

}
