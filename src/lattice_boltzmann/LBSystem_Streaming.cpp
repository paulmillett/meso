
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Execute the STREAMING step:
// -------------------------------------------------------------------------

void LBSystem::streamingStep()
{

   //	---------------------------------------
   //	Execute streaming algorithm:
   //	---------------------------------------

   for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            for (int fl=0; fl<nc; fl++) {
               for (int id=0; id<nn; id++) {
                  int inbr = nbrIndex(i,int(ex[id]),nx);
                  int jnbr = nbrIndex(j,int(ey[id]),ny);
                  f[inbr][jnbr][k][fl][id] = fprev[i][j][k][fl][id] +
                                             omega[i][j][k][fl][id];
               }
            }
         }
      }
   }

   //	---------------------------------------
   //	Fill in halo regions for the 'f' array:
   //	---------------------------------------

   // exchangeHaloDist();

}
