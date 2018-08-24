
# include "LBSystem.hpp"
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Parser for the initial condition:
// -------------------------------------------------------------------------

void LBSystem::parseInitialCondition()
{

	//	---------------------------------------
	//	Determine which solid IC to use:
	//	---------------------------------------

	string so_type = input_params("LBApp/initial_condition/solid","none");

	if (so_type == "uniform_channel") {initializeSolidChannel();}
	if (so_type == "input_file") {initializeSolidFromFile();}

	//	---------------------------------------
	//	Determine which fluid IC to use:
	//	---------------------------------------

	string ic_type = input_params("LBApp/initial_condition/type","binary_mixed");

	if (ic_type == "single_droplet") {initializeSingleDroplet();}
	if (ic_type == "binary_mixed") {initializeBinaryMixed();}
	if (ic_type == "binary_foam") {initializeBinaryFoam();}

}



// -------------------------------------------------------------------------
// Initialize SOLID CHANNEL: uniform width channel...
// -------------------------------------------------------------------------

void LBSystem::initializeSolidChannel()
{

	//	---------------------------------------
	//	set the channel walls:
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (j<3 || j>nx-3) {
					solid[i][j][k] = 1;
				}
			}
		}
	}
}



// -------------------------------------------------------------------------
// Initialize SOLID CHANNEL: read from input file...
// -------------------------------------------------------------------------

void LBSystem::initializeSolidFromFile()
{

	//	---------------------------------------
	//	open file:
	//	---------------------------------------

	ifstream infile;
	infile.open("solid.txt");

	//	---------------------------------------
	//	read data:
	//	---------------------------------------

	int nxin, nyin;
	infile >> nxin >> nyin;

	if (nxin == nxGlobal && nyin == ny) {
		for (int i=istr; i<iend+1; i++) {
	      for (int j=jstr; j<jend+1; j++) {
	         for (int k=kstr; k<kend+1; k++) {
					infile >> solid[i][j][k];
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Initialize BINARY MIXTURE: homogeneously mixed with some noise:
// -------------------------------------------------------------------------

void LBSystem::initializeBinaryMixed()
{

	//	---------------------------------------
	//	Set parameters:
	//	---------------------------------------

	double rhoA = input_params("LBApp/initial_condition/rhoA",0.5);
	double rhoB = input_params("LBApp/initial_condition/rhoB",0.5);

	//	---------------------------------------
   //	Set fluid density & velocity:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				double rhosum = 0.0;
				for (int fl=0; fl<nc; fl++) {
					double r1 = (double)rand()/RAND_MAX;
					if (fl == 0) rho[i][j][k][fl] = rhoA + 0.2*(r1 - 0.5);
					if (fl == 1) rho[i][j][k][fl] = rhoB + 0.2*(r1 - 0.5);
					rhosum += rho[i][j][k][fl];
					u[i][j][k][fl] = 0.0;
					v[i][j][k][fl] = 0.0;
				}
				rhoTot[i][j][k] = rhosum;
				uTot[i][j][k] = 0.0;
				vTot[i][j][k] = 0.0;
			}
		}
	}

	//	---------------------------------------
   //	Get the equilibrium distributions:
   //	---------------------------------------

	exchangeHaloRho();
	parseFluidForces();
	equilibriumDistribution();

	//	---------------------------------------
   //	Particle distributions assigned equil.
	// values:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				for (int fl=0; fl<nc; fl++) {
					for (int id=0; id<nn; id++) {
						f[i][j][k][fl][id] = feq[i][j][k][fl][id];
						fprev[i][j][k][fl][id] = f[i][j][k][fl][id];
					}
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Initialize a single droplet in a binary fluid:
// -------------------------------------------------------------------------

void LBSystem::initializeSingleDroplet()
{

	//	---------------------------------------
	//	Set parameters:
	//	---------------------------------------

	double xpos = input_params("LBApp/initial_condition/xpos",10.0);
	double ypos = input_params("LBApp/initial_condition/ypos",10.0);
	double rad = input_params("LBApp/initial_condition/rad",5.0);

	//	---------------------------------------
	//	assign a single droplet to system...
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
				double rx = abs(double(i) - xpos);
				double ry = abs(double(j) - ypos);
				if (rx > 0.5*nx) rx = nx - rx;
				if (ry > 0.5*ny) ry = ny - ry;
				double r = sqrt(rx*rx + ry*ry);
				if (r <= rad) {
					rho[i][j][k][0] = 1.15;
					rho[i][j][k][1] = 0.05;
				}
				else {
					rho[i][j][k][0] = 0.05;
					rho[i][j][k][1] = 1.15;
				}
				rhoTot[i][j][k] = rho[i][j][k][0] + rho[i][j][k][1];
				uTot[i][j][k] = 0.0;
				vTot[i][j][k] = 0.0;
			}
		}
	}

	//	---------------------------------------
   //	Get the equilibrium distributions:
   //	---------------------------------------

	exchangeHaloRho();
	parseFluidForces();
	equilibriumDistribution();

	//	---------------------------------------
   //	Particle distributions assigned equil.
	// values:
   //	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				for (int fl=0; fl<nc; fl++) {
					for (int id=0; id<nn; id++) {
						f[i][j][k][fl][id] = feq[i][j][k][fl][id];
						fprev[i][j][k][fl][id] = f[i][j][k][fl][id];
					}
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Initialize BINARY FOAM: a Voronoi tesselation is used for a dense droplet
// array:
// -------------------------------------------------------------------------

void LBSystem::initializeBinaryFoam()
{

	//	---------------------------------------
	//	Set parameters:
	//	---------------------------------------

	int ndrops = input_params("LBApp/initial_condition/ndrops",10);

	double dcx[50];  // arrays that store x- and y-positions of
	double dcy[50];  // droplets

	int** flag = new int*[nx+1];  // these are arrays
	int**  tag = new int*[nx+1];  // that define the Voronoi
	for (int i=0; i<nx+1; i++) {  // tesselation
		flag[i] = new int[ny+1];
		 tag[i] = new int[ny+1];
	}

	//	---------------------------------------
	//	FIRST: define droplet centers...
	//	---------------------------------------

	for (int i=1; i<ndrops+1; i++) {
		droplet_attempt:
		double r1 = (double)rand()/RAND_MAX;
		double r2 = (double)rand()/RAND_MAX;
		dcx[i] = r1*nx + 1.0;
		dcy[i] = r2*ny + 1.0;
		for (int j=1; j<i; j++) {
			double rx = abs(dcx[i] - dcx[j]);
			double ry = abs(dcy[i] - dcy[j]);
			if (rx > 0.5*nx) rx = nx - rx;
			if (ry > 0.5*ny) ry = ny - ry;
			double r = sqrt(rx*rx + ry*ry);
			if (r < 10.0) goto droplet_attempt; // go back to beginning of loop;
		}
	}

	//	---------------------------------------
	//	SECOND: flag grid w/ voronoi tesselation
	//	---------------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			double r1st = 2.0*nx;
			int idrop = 0;
			for (int k=1; k<ndrops+1; k++) {
				double rx = abs(double(i) - dcx[k]);
				double ry = abs(double(j) - dcy[k]);
				if (rx > 0.5*nx) rx = nx - rx;
				if (ry > 0.5*ny) ry = ny - ry;
				double r = sqrt(rx*rx + ry*ry);
				if (r < r1st) {
					r1st = r;
					idrop = k;
				}
			}
			flag[i][j] = idrop;
		}
	}

	//	---------------------------------------
	//	THIRD: iteratively 'expand' boundaries
	// between drops using 'tag' and 'flag'
	// arrays...
	//	---------------------------------------

	for (int k=1; k<3; k++) {
		// set tag = flag:
		for (int i=1; i<nx+1; i++) {
			for (int j=1; j<ny+1; j++) {
				tag[i][j] = flag[i][j];
			}
		}
		// tag a node if it's nabor is not the same droplet:
		for (int i=1; i<nx+1; i++) {
			for (int j=1; j<ny+1; j++) {
				for (int id=1; id<5; id++) {
					int inbr = nbrIndex(i,int(ex[id]),nx);
					int jnbr = nbrIndex(j,int(ey[id]),ny);
					if (flag[inbr][jnbr] != flag[i][j] || solid[i][j][k] == 1) {
						tag[i][j] = 0;
						break;
					}
				}
			}
		}
		// set flag = tag:
		for (int i=1; i<nx+1; i++) {
			for (int j=1; j<ny+1; j++) {
				flag[i][j] = tag[i][j];
			}
		}
	}

	//	---------------------------------------
	//	FOURTH: assign fluid densities...
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				if (solid[i][j][k] == 1) continue;
				if (flag[i][j] > 0) {
					rho[i][j][k][0] = 1.15;
					rho[i][j][k][1] = 0.05;
				}
				else {
					rho[i][j][k][0] = 0.05;
					rho[i][j][k][1] = 1.15;
				}
				for (int fl=0; fl<nc; fl++) {
					u[i][j][k][fl] = 0.0;
					v[i][j][k][fl] = 0.0;
				}
				rhoTot[i][j][k] = rho[i][j][k][0] + rho[i][j][k][1];
				uTot[i][j][k] = 0.0;
				vTot[i][j][k] = 0.0;
			}
		}
	}

	//	---------------------------------------
	//	De-allocate arrays:
	//	---------------------------------------

	for (int i=0; i<nx+1; i++) {
		delete[] tag[i];
		delete[] flag[i];
	}
	delete[] tag;
	delete[] flag;

	//	---------------------------------------
	//	Get the equilibrium distributions:
	//	---------------------------------------

	parseFluidForces();
	equilibriumDistribution();

	//	---------------------------------------
	//	Particle distributions assigned equil.
	// values:
	//	---------------------------------------

	for (int i=istr; i<iend+1; i++) {
      for (int j=jstr; j<jend+1; j++) {
         for (int k=kstr; k<kend+1; k++) {
				for (int fl=0; fl<nc; fl++) {
					for (int id=0; id<nn; id++) {
						f[i][j][k][fl][id] = feq[i][j][k][fl][id];
						fprev[i][j][k][fl][id] = f[i][j][k][fl][id];
					}
				}
			}
		}
	}

}
