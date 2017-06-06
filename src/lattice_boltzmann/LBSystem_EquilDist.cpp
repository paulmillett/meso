
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Calculate the equilibrium particle distributions for each component:
// NOTE: the below algorithm uses the BGK approximation.
// -------------------------------------------------------------------------

void LBSystem::equilibriumDistribution()
{
   for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            if (solid[i][j][k] == 1) continue;
            for (int fl=0; fl<nc; fl++) {
               // modified velocities due to forces...
               double u_F = uTot[i][j][k]+ tau*fx[i][j][k][fl]/rho[i][j][k][fl];
               double v_F = vTot[i][j][k]+ tau*fy[i][j][k][fl]/rho[i][j][k][fl];
               double uv2 = u_F*u_F + v_F*v_F;
               for (int id=0; id<nn; id++) {
                  double evel = ex[id]*u_F + ey[id]*v_F;
                  feq[i][j][k][fl][id] = rho[i][j][k][fl]*wa[id]*(1.0 +
                                         3.0*evel + 4.5*evel*evel - 1.5*uv2);
               }
            }
         }
      }
   }
}
