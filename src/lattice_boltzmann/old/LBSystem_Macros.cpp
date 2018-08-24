
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Macros: calculate total velocity and density:
// -------------------------------------------------------------------------

void LBSystem::calculateMacros()
{

	//	---------------------------------------
	//	Update fprev = f:
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            for (int fl=0; fl<nc; fl++) {
               for (int id=0; id<nn; id++) {
						fprev[i][j][k][fl][id] = f[i][j][k][fl][id];
					}
				}
			}
		}
	}

   //	---------------------------------------
   //	Calculate local density & velocity:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            if (solid[i][j][k] == 1) continue;
            double sumrho = 0.0;
            double sumrhou = 0.0;
            double sumrhov = 0.0;
            for (int fl=0; fl<nc; fl++) {
               double sum = 0.0;
               double sumx = 0.0;
               double sumy = 0.0;
               for (int id=0; id<nn; id++) {
                  sum += f[i][j][k][fl][id];
                  sumx += f[i][j][k][fl][id]*ex[id];
                  sumy += f[i][j][k][fl][id]*ey[id];
               }
               rho[i][j][k][fl] = sum;
               u[i][j][k][fl] = sumx/rho[i][j][k][fl];
               v[i][j][k][fl] = sumy/rho[i][j][k][fl];
               sumrho += rho[i][j][k][fl];
               sumrhou += rho[i][j][k][fl]*u[i][j][k][fl];
               sumrhov += rho[i][j][k][fl]*v[i][j][k][fl];
            }
            rhoTot[i][j][k] = sumrho;
            uTot[i][j][k] = sumrhou/rhoTot[i][j][k];
            vTot[i][j][k] = sumrhov/rhoTot[i][j][k];
         }
      }
   }

	//	---------------------------------------
   //	Fill in halo regions for the 'rho' array:
   //	---------------------------------------

	// exchangeHaloRho();

}
