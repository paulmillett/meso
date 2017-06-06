
# include "LBSystem.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Parser for chemical potentials:
// -------------------------------------------------------------------------

void LBSystem::parseFluidForces()
{

	//	---------------------------------------
	//	Determine which fluid-force model to use:
	//	---------------------------------------

	string force_type = input_params("LBApp/fluid_force/type","mcmp");

	if (force_type == "mcmp") {calculateFluidForcesMCMP();}
   if (force_type == "mcmp_foam") {calculateFluidForcesMCMPfoam();}

	//	---------------------------------------
	//	Determine which body-force model to use:
	//	---------------------------------------

	string bodyf_type = input_params("LBApp/body_force/type","none");

	if (bodyf_type == "uniform_x_dir") {addBodyForceXDir();}

}



// -------------------------------------------------------------------------
// Calculate fluid forces for each component using:
//	                          Multi-component, Multi-phase model of Shan-Chen
// -------------------------------------------------------------------------

void LBSystem::calculateFluidForcesMCMP()
{

	//	---------------------------------------
	//	Get needed parameters:
	//	---------------------------------------

	double G = input_params("LBApp/fluid_force/G",3.5);  // fluid-fluid interaction param.
	double Gs[2];
	Gs[0] = 2.0*G;
	Gs[1] = 0.3*G;

   //	---------------------------------------
   //	Zero out forces:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            for (int fl=0; fl<nc; fl++) {
					fx[i][j][k][fl] = 0.0;
               fy[i][j][k][fl] = 0.0;
				}
			}
		}
	}

   //	---------------------------------------
   //	calculate intermolecular forces:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
            for (int fl=0; fl<nc; fl++) {
               int op = 1 - fl;  // opposite fluid index
               double Gsumx = 0.0;
               double Gsumy = 0.0;
					double Ssumx = 0.0;
					double Ssumy = 0.0;
               for (int id=0; id<nn; id++) {
						int inbr = nbrIndex(i,int(ex[id]),nx);
                  int jnbr = nbrIndex(j,int(ey[id]),ny);
                  Gsumx += wa[id]*psi(rho[inbr][jnbr][k][op])*ex[id];  // fluid-fluid
                  Gsumy += wa[id]*psi(rho[inbr][jnbr][k][op])*ey[id];
						Ssumx += wa[id]*psi(solid[inbr][jnbr][k])*ex[id];    // solid-fluid
						Ssumy += wa[id]*psi(solid[inbr][jnbr][k])*ey[id];
               }
               fx[i][j][k][fl] = -G*psi(rho[i][j][k][fl])*Gsumx - Gs[fl]*psi(rho[i][j][k][fl])*Ssumx;
               fy[i][j][k][fl] = -G*psi(rho[i][j][k][fl])*Gsumy - Gs[fl]*psi(rho[i][j][k][fl])*Ssumy;
            }
         }
      }
   }

}



// -------------------------------------------------------------------------
// Calculate fluid forces for each component using:
//	                          Multi-component, Multi-phase model of Sbragaglia
// -------------------------------------------------------------------------

void LBSystem::calculateFluidForcesMCMPfoam()
{

	//	---------------------------------------
	//	Get needed parameters:
	//	---------------------------------------

	double gAB = input_params("LBApp/fluid_force/gAB",0.8265);  // A-B inter. forces
	double gAAshort = input_params("LBApp/fluid_force/gAAshort",-15.0); // A-A short-range
	double gBBshort = input_params("LBApp/fluid_force/gBBshort",-14.0); // B-B short-range
	double gAAlong = input_params("LBApp/fluid_force/gAAlong",14.1); // A-A long-range
	double gBBlong = input_params("LBApp/fluid_force/gBBlong",13.1); // B-B long-range

	double gSH[2]; gSH[0] = gAAshort; gSH[1] = gBBshort;  // short-range force consts.
	double gLO[2]; gLO[0] = gAAlong;  gLO[1] = gBBlong;   // long-range  force consts.
	double gFS[2]; gFS[0] = 1.8*gAB;  gFS[1] = 0.5*gAB;   // fluid-solid force consts.

   //	---------------------------------------
   //	Zero out forces:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            for (int fl=0; fl<nc; fl++) {
					fx[i][j][k][fl] = 0.0;
               fy[i][j][k][fl] = 0.0;
				}
			}
		}
	}

   //	---------------------------------------
   //	calculate A-B intermolecular forces:
	// {these will be repulsive}
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
            for (int fl=0; fl<nc; fl++) {
               int op = 1 - fl;  // opposite fluid index
               double Gsumx = 0.0;
               double Gsumy = 0.0;
					double Ssumx = 0.0;
					double Ssumy = 0.0;
               for (int id=1; id<9; id++) {
						int inbr = nbrIndex(i,int(ex[id]),nx);
                  int jnbr = nbrIndex(j,int(ey[id]),ny);
                  Gsumx += wa[id]*rho[inbr][jnbr][k][op]*ex[id];  // fluid-fluid
                  Gsumy += wa[id]*rho[inbr][jnbr][k][op]*ey[id];
						Ssumx += wa[id]*solid[inbr][jnbr][k]*ex[id];    // fluid-solid
						Ssumy += wa[id]*solid[inbr][jnbr][k]*ey[id];
               }
               fx[i][j][k][fl] = -gAB*rho[i][j][k][fl]*Gsumx - gFS[fl]*rho[i][j][k][fl]*Ssumx;
               fy[i][j][k][fl] = -gAB*rho[i][j][k][fl]*Gsumy - gFS[fl]*rho[i][j][k][fl]*Ssumy;
            }
         }
      }
   }

	//	---------------------------------------
   //	calculate short-range A-A & B-B forces:
	// {these will be attractive}
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
				for (int fl=0; fl<nc; fl++) {
					double Gsumx = 0.0;
					double Gsumy = 0.0;
					for (int id=1; id<9; id++) {
						int inbr = nbrIndex(i,int(ex[id]),nx);
                  int jnbr = nbrIndex(j,int(ey[id]),ny);
						Gsumx += wa[id]*psi2(rho[inbr][jnbr][k][fl])*ex[id];
						Gsumy += wa[id]*psi2(rho[inbr][jnbr][k][fl])*ey[id];
					}
					fx[i][j][k][fl] -= gSH[fl]*psi2(rho[i][j][k][fl])*Gsumx;
					fy[i][j][k][fl] -= gSH[fl]*psi2(rho[i][j][k][fl])*Gsumy;
				}
			}
		}
	}

	//	---------------------------------------
   //	calculate long-range A-A & B-B forces:
	// {these will be repulsive}
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
				for (int fl=0; fl<nc; fl++) {
					double Gsumx = 0.0;
					double Gsumy = 0.0;
					for (int id=1; id<25; id++) {
						int inbr = nbrIndex(i,int(ex[id]),nx);
                  int jnbr = nbrIndex(j,int(ey[id]),ny);
						Gsumx += pa[id]*psi2(rho[inbr][jnbr][k][fl])*ex[id];
						Gsumy += pa[id]*psi2(rho[inbr][jnbr][k][fl])*ey[id];
					}
					fx[i][j][k][fl] -= gLO[fl]*psi2(rho[i][j][k][fl])*Gsumx;
					fy[i][j][k][fl] -= gLO[fl]*psi2(rho[i][j][k][fl])*Gsumy;
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Potential function of the fluid density:
// -------------------------------------------------------------------------

double LBSystem::psi(double rho_in)
{
   return (1.0 - exp(-rho_in));
}



// -------------------------------------------------------------------------
// Alternative potential function of the fluid density:
// -------------------------------------------------------------------------

double LBSystem::psi2(double rho_in)
{
	double rho0 = 0.7;
	return rho0*(1.0 - exp(-rho_in/rho0));
}



// -------------------------------------------------------------------------
// neighbor index calculation:
// -------------------------------------------------------------------------

int LBSystem::nbrIndex(int i, int di, int ni)
{
   int val = i + di;
   if (val < 1)  val = val + ni;
   if (val > ni) val = val - ni;
   return val;
}



// -------------------------------------------------------------------------
// Add body force that is applied in the x-direction:
// -------------------------------------------------------------------------

void LBSystem::addBodyForceXDir()
{

	//	---------------------------------------
	//	Get needed parameters:
	//	---------------------------------------

	double xBF = input_params("LBApp/body_force/xBF",0.0);

	//	---------------------------------------
	//	add force:
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
				for (int fl=0; fl<nc; fl++) {
					fx[i][j][k][fl] += xBF;
				}
			}
		}
	}
}
