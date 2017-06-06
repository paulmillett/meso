
# include "ZoneAnnealing.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

ZoneAnnealing::ZoneAnnealing(const Params& p,
                             const vector<FDVector>& cin) : c(cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   nc = p.nc;
   vzone = p.input_params("CHApp/mobility/vzone",1.0);
   wzone = p.input_params("CHApp/mobility/wzone",20.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

ZoneAnnealing::~ZoneAnnealing()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

double ZoneAnnealing::mobFunc(int ndx, int n, int i, int j, int k)
{
   return 1.0;  // need to define later
}
