
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Execute the COLLISION step:
// -------------------------------------------------------------------------

void LBSystem::collisionStep()
{
   for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            if (solid[i][j][k] == 1) continue;
            for (int fl=0; fl<nc; fl++) {
               for (int id=0; id<nn; id++) {
                  double fdiff = fprev[i][j][k][fl][id] - feq[i][j][k][fl][id];
                  omega[i][j][k][fl][id] = -fdiff/tau;
               }
            }
         }
      }
   }
}
