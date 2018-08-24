
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Exchange halo data for the 'rho' array:
// -------------------------------------------------------------------------

void LBSystem::exchangeHaloRho()
{

   //	---------------------------------------
   //	exchange in x-direction:
   //	---------------------------------------

   for (int j=0; j<ny+halo+1; j++) {
      for (int k=1; k<nz+1; k++) {
         for (int fl=0; fl<nc; fl++) {
            rho[0][j][k][fl] = rho[nx][j][k][fl];
            rho[nx+1][j][k][fl] = rho[1][j][k][fl];
         }
      }
   }

   //	---------------------------------------
   //	exchange in y-direction:
   //	---------------------------------------

   for (int i=0; i<nx+halo+1; i++) {
      for (int k=1; k<nz+1; k++) {
         for (int fl=0; fl<nc; fl++) {
            rho[i][0][k][fl] = rho[i][ny][k][fl];
            rho[i][ny+1][k][fl] = rho[i][1][k][fl];
         }
      }
   }

}



// -------------------------------------------------------------------------
// Exchange halo data for the 'f' array:
// -------------------------------------------------------------------------

void LBSystem::exchangeHaloDist()
{

   //	---------------------------------------
   //	exchange in x-direction:
   //	---------------------------------------

   // for (int j=0; j<ny+halo+1; j++) {
   //    for (int k=1; k<nz+1; k++) {
   //       for (int fl=0; fl<nc; fl++) {
   //          for (int id=0; id<nn; id++) {
   //             f[nx][j][k][fl][id] = f[0][j][k][fl][id];
   //             f[1][j][k][fl][id] = f[nx+1][j][k][fl][id];
   //          }
   //       }
   //    }
   // }

   // for (int j=jstr; j<jend+1; j++) {
   //    for (int k=1; k<nz+1; k++) {
   //       for (int fl=0; fl<nc; fl++) {
   //          for (int id=0; id<nn; id++) {
   //             // update particles streaming to the right...
   //             int isrce = iend + 1;
   //             int idest = istr;
   //             int jdest = j; //nbrIndex(j,int(ey[id]),ny);
   //             if (ex[id] > 0.0) {
   //                f[idest][jdest][fl][id] = f[isrce][j][fl][id];
   //                if (idest == 1 && jdest == 1 && fl==0 && id==5) cout << f[idest][jdest][k][fl][id] << " exch" << endl;
   //             }
   //
   //
   //
   //
   //             // update particles streaming to the left...
   //             isrce = istr - 1;
   //             idest = iend;
   //             jdest = j; //nbrIndex(j,int(ey[id]),ny);
   //             if (ex[id] < 0.0) f[idest][jdest][fl][id] = f[isrce][j][fl][id];
   //          }
   //       }
   //    }
   // }

   //	---------------------------------------
   //	exchange in y-direction:
   //	---------------------------------------

   // for (int i=0; i<nx+halo+1; i++) {
   //    for (int k=1; k<nz+1; k++) {
   //       for (int fl=0; fl<nc; fl++) {
   //          for (int id=0; id<nn; id++) {
   //             f[i][ny][k][fl][id] = f[i][0][k][fl][id];
   //             f[i][1][k][fl][id] = f[i][ny+1][k][fl][id];
   //          }
   //       }
   //    }
   // }

   // for (int i=istr; i<iend+1; i++) {
   //    for (int k=1; k<nz+1; k++) {
   //       for (int fl=0; fl<nc; fl++) {
   //          for (int id=0; id<nn; id++) {
   //             // update particles streaming to the up...
   //             int jsrce = jend + 1;
   //             int jdest = jstr;
   //             int idest = i; //nbrIndex(i,int(ex[id]),nx);
   //             if (ey[id] > 0.0) f[idest][jdest][fl][id] = f[i][jsrce][fl][id];
   //             // update particles streaming to the down...
   //             jsrce = jstr - 1;
   //             jdest = jend;
   //             idest = i; //nbrIndex(i,int(ex[id]),nx);
   //             if (ey[id] < 0.0) f[idest][jdest][fl][id] = f[i][jsrce][fl][id];
   //          }
   //       }
   //    }
   // }

}
