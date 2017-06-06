
# include "DiblockCopolymer.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

DiblockCopolymer::DiblockCopolymer(const Params& p,
                                   const vector<FDVector>& cin) : c(cin)
{

   // ----------------------------------------
   //	Establish the parameters:
   // ----------------------------------------

   nc = p.nc;
   A0 = p.input_params("CHApp/source/A0",1.0);
   c0_init_mean = p.input_params("CHApp/initial_condition/c0_init_mean",0.5);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

DiblockCopolymer::~DiblockCopolymer()
{

}



// -------------------------------------------------------------------------
// Function to calculate mu:
// -------------------------------------------------------------------------

double DiblockCopolymer::srcFunc(int ndx, int n, int i, int j, int k)
{
   double cc = c[n].getValue(ndx);
   return -A0*(cc - c0_init_mean);
}
