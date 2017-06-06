
# include "LBSystem.hpp"
# include <string>
# include <iomanip>
# include <sstream>
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

LBSystem::LBSystem(const GetPot& in_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

	input_params = in_params;
	nc = input_params("LBApp/nc",1);
	nn = input_params("LBApp/nn",9);
   nxGlobal = input_params("Domain/nx",1);
   ny = input_params("Domain/ny",1);
   nz = input_params("Domain/nz",1);
   dx = input_params("Domain/dx",1.0);
   dy = input_params("Domain/dy",1.0);
   dz = input_params("Domain/dz",1.0);
   dt = input_params("Time/dt",1.0);

   nx = nxGlobal;   // only for non-parallel jobs - testing
	xOff = 0;        // "  "

	//	---------------------------------------
   //	Get MPI rank:
   //	---------------------------------------

   rank = MPI::COMM_WORLD.Get_rank();

	//	---------------------------------------
   //	Initialize parameters:
   //	---------------------------------------

	tau = 1.0;
	mu = (2.0*tau - 1.0)*dx*dx/(6.0*dt);

	//	---------------------------------------
	//	Define lattice information:
	//	---------------------------------------

	lattice_type = input_params("LBApp/lattice","D2Q9");

	// D2Q9 lattice type:

	if (lattice_type == "D2Q9") {
		nn = 9;
		nnn = 9;
		halo = 1;
		ex = new double[nn];
		ey = new double[nn];
		wa = new double[nn];

		// lattice vectors...
		ex[0]=0; ex[1]=1; ex[2]=0; ex[3]=-1; ex[4]=0; ex[5]=1; ex[6]=-1; ex[7]=-1; ex[8]=1;
		ey[0]=0; ey[1]=0; ey[2]=1; ey[3]=0; ey[4]=-1; ey[5]=1; ey[6]=1; ey[7]=-1; ey[8]=-1;

		// weight factors...
		wa[0]=4./9.; wa[1]=1./9.; wa[2]=1./9.; wa[3]=1./9.; wa[4]=1./9.;
		wa[5]=1./36.; wa[6]=1./36.; wa[7]=1./36.; wa[8]=1./36.;

	}

	// D2Q25 lattice type:

	if (lattice_type == "D2Q25") {
		nn = 9;
		nnn = 25;
		halo = 1;  //2;
		ex = new double[25];
		ey = new double[25];
		wa = new double[9];
		pa = new double[25];

		// lattice vectors...
		ex[0]=0; ex[1]=1; ex[2]=0; ex[3]=-1; ex[4]=0; ex[5]=1; ex[6]=-1; ex[7]=-1; ex[8]=1;
		ey[0]=0; ey[1]=0; ey[2]=1; ey[3]=0; ey[4]=-1; ey[5]=1; ey[6]=1; ey[7]=-1; ey[8]=-1;

		ex[9]=2; ex[10]=0; ex[11]=-2; ex[12]=0;
		ey[9]=0; ey[10]=2; ey[11]=0; ey[12]=-2;

		ex[13]=2; ex[14]=1; ex[15]=-1; ex[16]=-2; ex[17]=-2; ex[18]=-1; ex[19]= 1; ex[20]= 2;
		ey[13]=1; ey[14]=2; ey[15]= 2; ey[16]= 1; ey[17]=-1; ey[18]=-2; ey[19]=-2; ey[20]=-1;

		ex[21]=2; ex[22]=-2; ex[23]=-2; ex[24]= 2;
		ey[21]=2; ey[22]= 2; ey[23]=-2; ey[24]=-2;

		// weight factors for 1st shell neighbors...
		wa[0]=4./9.; wa[1]=1./9.; wa[2]=1./9.; wa[3]=1./9.; wa[4]=1./9.;
		wa[5]=1./36.; wa[6]=1./36.; wa[7]=1./36.; wa[8]=1./36.;

		// weight factors for 2nd shell neighbors...
		double p0  = 247.0/420.0;
		double p1  = 4.0/63.0;
		double p5  = 4.0/135.0;
		double p9  = 1.0/180.0;
		double p13 = 2.0/945.0;
		double p21 = 1.0/15120.0;

		pa[0] = p0; pa[1] = p1; pa[2] = p1; pa[3] = p1; pa[4] = p1;
		pa[5] = p5; pa[6] = p5; pa[7] = p5; pa[8] = p5;
		pa[9] = p9; pa[10] = p9; pa[11] = p9; pa[12] = p9;
		pa[13] = p13; pa[14] = p13; pa[15] = p13; pa[16] = p13;
		pa[17] = p13; pa[18] = p13; pa[19] = p13; pa[20] = p13;
		pa[21] = p21; pa[22] = p21; pa[23] = p21; pa[24] = p21;

	}

	//	---------------------------------------
	//	Define start & end indices for data:
	//	---------------------------------------

	istr = 0 + halo;
	iend = istr + nx - 1;
	jstr = 0 + halo;
	jend = jstr + ny - 1;
	kstr = 1;
	kend = 1;

   //	---------------------------------------
   //	Allocate data arrays:
	// NOTE: array sizes are expanded due to halo
   //	---------------------------------------

   solid = new int** [nx+2*halo];
	nbr = new int**** [nx+2*halo];
   rhoTot = new double** [nx+2*halo];
   uTot = new double** [nx+2*halo];
   vTot = new double** [nx+2*halo];
   rho = new double*** [nx+2*halo];
   u = new double*** [nx+2*halo];
   v = new double*** [nx+2*halo];
   fx = new double*** [nx+2*halo];
   fy = new double*** [nx+2*halo];
   omega = new double**** [nx+2*halo];
   f = new double**** [nx+2*halo];
   feq = new double**** [nx+2*halo];
   fprev = new double**** [nx+2*halo];
	for (int i=0; i<nx+2; i++) {
		solid[i] = new int* [ny+2*halo];
		nbr[i] = new int*** [ny+2*halo];
      rhoTot[i] = new double* [ny+2*halo];
      uTot[i] = new double* [ny+2*halo];
      vTot[i] = new double* [ny+2*halo];
      rho[i] = new double** [ny+2*halo];
      u[i] = new double** [ny+2*halo];
      v[i] = new double** [ny+2*halo];
      fx[i] = new double** [ny+2*halo];
      fy[i] = new double** [ny+2*halo];
      omega[i] = new double*** [ny+2*halo];
      f[i] = new double*** [ny+2*halo];
      feq[i] = new double*** [ny+2*halo];
      fprev[i] = new double*** [ny+2*halo];
      for (int j=0; j<ny+2; j++) {
         solid[i][j] = new int[nz+2];
			nbr[i][j] = new int**[nz+2];
         rhoTot[i][j] = new double[nz+2];
         uTot[i][j] = new double[nz+2];
         vTot[i][j] = new double[nz+2];
         rho[i][j] = new double* [nz+2];
         u[i][j] = new double* [nz+2];
         v[i][j] = new double* [nz+2];
         fx[i][j] = new double* [nz+2];
         fy[i][j] = new double* [nz+2];
         omega[i][j] = new double** [nz+2];
         f[i][j] = new double** [nz+2];
         feq[i][j] = new double** [nz+2];
         fprev[i][j] = new double** [nz+2];
         for (int k=0; k<nz+2; k++) {
				nbr[i][j][k] = new int*[nnn];
				rho[i][j][k] = new double[nc];
            u[i][j][k] = new double[nc];
            v[i][j][k] = new double[nc];
            fx[i][j][k] = new double[nc];
            fy[i][j][k] = new double[nc];
            omega[i][j][k] = new double* [nc];
            f[i][j][k] = new double* [nc];
            feq[i][j][k] = new double* [nc];
            fprev[i][j][k] = new double* [nc];
            for (int m=0; m<nc; m++) {
               omega[i][j][k][m] = new double [nn];
               f[i][j][k][m] = new double [nn];
               feq[i][j][k][m] = new double [nn];
               fprev[i][j][k][m] = new double [nn];
            }
				for (int m=0; m<nnn; m++) {
					nbr[i][j][k][m] = new int [3];
				}
         }
      }
	}

	//	---------------------------------------
	//	Define neighbor indices for each point
	// in each direction:
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
            for (int id=0; id<nnn; id++) {
					nbr[i][j][k][id][0] = nbrIndex(i,int(ex[id]),nx);
					nbr[i][j][k][id][1] = nbrIndex(j,int(ey[id]),ny);
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBSystem::~LBSystem()
{

   //	---------------------------------------
   //	Deallocate data arrays:
   //	---------------------------------------

   for (int i=0; i<nx+2; i++) {
      for (int j=0; j<ny+2; j++) {
         delete[] solid[i][j];
         delete[] rhoTot[i][j];
         delete[] uTot[i][j];
         delete[] vTot[i][j];
         for (int k=0; k<nz+2; k++) {
            delete[] rho[i][j][k];
            delete[] u[i][j][k];
            delete[] v[i][j][k];
            delete[] fx[i][j][k];
            delete[] fy[i][j][k];
            for (int m=0; m<nc; m++) {
               delete[] omega[i][j][k][m];
               delete[] f[i][j][k][m];
               delete[] feq[i][j][k][m];
               delete[] fprev[i][j][k][m];
            }
				for (int m=0; m<nnn; m++) {
					delete[] nbr[i][j][k][m];
				}
				delete[] nbr[i][j][k];
            delete[] omega[i][j][k];
            delete[] f[i][j][k];
            delete[] feq[i][j][k];
            delete[] fprev[i][j][k];
         }
			delete[] nbr[i][j];
			delete[] rho[i][j];
         delete[] u[i][j];
         delete[] v[i][j];
         delete[] fx[i][j];
         delete[] fy[i][j];
         delete[] omega[i][j];
         delete[] f[i][j];
         delete[] feq[i][j];
         delete[] fprev[i][j];
      }
		delete[] nbr[i];
		delete[] solid[i];
      delete[] rhoTot[i];
      delete[] uTot[i];
      delete[] vTot[i];
      delete[] rho[i];
      delete[] u[i];
      delete[] v[i];
      delete[] fx[i];
      delete[] fy[i];
      delete[] omega[i];
      delete[] f[i];
      delete[] feq[i];
      delete[] fprev[i];
	}
	delete[] nbr;
	delete[] solid;
   delete[] rhoTot;
   delete[] uTot;
   delete[] vTot;
   delete[] rho;
   delete[] u;
   delete[] v;
   delete[] fx;
   delete[] fy;
   delete[] omega;
   delete[] f;
   delete[] feq;
   delete[] fprev;

	delete[] ex;
	delete[] ey;
	delete[] wa;

}
