
# include "SimpleDoubleWell.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

SimpleDoubleWell::SimpleDoubleWell(const Params& p,
                                   const vector<FDVector>& cin) : c(cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   nc = p.nc;
   w = p.input_params("CHApp/free_energy/w",1.0);
   kap = p.input_params("CHApp/free_energy/kap",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

SimpleDoubleWell::~SimpleDoubleWell()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

double SimpleDoubleWell::muFunc(int ndx, int n, int i, int j, int k)
{
   double cc = c[n].getValue(ndx);
   double df = w*(4*cc*cc*cc - 6*cc*cc + 2*cc);
   double lap = c[n].calculateLaplacian(ndx);
   return df - kap*lap;
}
