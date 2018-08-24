
# include "LBSystem.hpp"
using namespace std;


// -------------------------------------------------------------------------
// Bounce-back conditions for fluid interaction with solid walls:
// -------------------------------------------------------------------------

void LBSystem::bounceBack()
{
	double temp;
	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            if (solid[i][j][k] == 1) {
               for (int fl=0; fl<nc; fl++) {
                  temp = f[i][j][k][fl][1]; f[i][j][k][fl][1] = f[i][j][k][fl][3]; f[i][j][k][fl][3] = temp;
                  temp = f[i][j][k][fl][2]; f[i][j][k][fl][2] = f[i][j][k][fl][4]; f[i][j][k][fl][4] = temp;
                  temp = f[i][j][k][fl][5]; f[i][j][k][fl][5] = f[i][j][k][fl][7]; f[i][j][k][fl][7] = temp;
                  temp = f[i][j][k][fl][6]; f[i][j][k][fl][6] = f[i][j][k][fl][8]; f[i][j][k][fl][8] = temp;
               }
            }
         }
      }
   }
}
