
# include "mcmp.hpp"
# include <fstream>
# include <iostream>
# include <iomanip>
# include <sstream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

mcmp::mcmp(const GetPot& in_params)
{

   //	---------------------------------------
   //	Assign variables from 'input_params':
   //	---------------------------------------

   input_params = in_params;
   nc = input_params("LBApp/nc",2);
   nn = input_params("LBApp/nn",9);
   NX = input_params("Domain/nx",1);
   NY = input_params("Domain/ny",1);
   NZ = input_params("Domain/nz",1);
   dx = input_params("Domain/dx",1.0);
   dy = input_params("Domain/dy",1.0);
   dz = input_params("Domain/dz",1.0);
   dt = input_params("Time/dt",1.0);
   G  = input_params("LBApp/fluid_force/G",3.5);  // fluid-fluid interaction param.

   nx = NX;   // only for non-parallel jobs - testing
   ny = NY;
   nz = NZ;

   deli = (nz+2)*(ny+2);
   delj = (nz+2);
   delk = 1;

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
      for (int i=0; i<nn; i++) {
   		ex.push_back(0.0);
         ey.push_back(0.0);
         wa.push_back(0.0);
   	}

      // lattice vectors...
      ex[0]=0; ex[1]=1; ex[2]=0; ex[3]=-1; ex[4]=0; ex[5]=1; ex[6]=-1; ex[7]=-1; ex[8]=1;
      ey[0]=0; ey[1]=0; ey[2]=1; ey[3]=0; ey[4]=-1; ey[5]=1; ey[6]=1; ey[7]=-1; ey[8]=-1;

      // weight factors...
      wa[0]=4./9.; wa[1]=1./9.; wa[2]=1./9.; wa[3]=1./9.; wa[4]=1./9.;
      wa[5]=1./36.; wa[6]=1./36.; wa[7]=1./36.; wa[8]=1./36.;

   }

   //	---------------------------------------
   //	Allocate vector sizes:
   //	---------------------------------------

   nxyz = (nx+2)*(ny+2)*(nz+2);  // # of grid points (including halo)
   nxyzc = nxyz*nc;
   nxyzcn = nxyzc*nn;

   for (int i=0; i<nxyz; i++) {
      rhoTot.push_back(0.0);
      uTot.push_back(0.0);
      vTot.push_back(0.0);
   }

   for (int i=0; i<nxyzc; i++) {
      rho.push_back(0.0);
      u.push_back(0.0);
      v.push_back(0.0);
      fx.push_back(0.0);
      fy.push_back(0.0);
   }

   for (int i=0; i<nxyzcn; i++) {
      omega.push_back(0.0);
      f.push_back(0.0);
      feq.push_back(0.0);
      fprev.push_back(0.0);
   }

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

mcmp::~mcmp()
{

}



// -------------------------------------------------------------------------
// Initialize the Lattice-Boltzmann system:
// -------------------------------------------------------------------------

void mcmp::initializeMCMP()
{

   //	---------------------------------------
	//	Set parameters:
	//	---------------------------------------

	double rhoA = input_params("LBApp/initial_condition/rhoA",0.5);
	double rhoB = input_params("LBApp/initial_condition/rhoB",0.5);

	//	---------------------------------------
   //	Set fluid density & velocity:
   //	---------------------------------------

	for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {
            int ndx = i*deli + j*delj + k*delk;
				double rhosum = 0.0;
				for (int c=0; c<nc; c++) {
               int ndxc = ndx*nc + c;
					double r1 = (double)rand()/RAND_MAX;
					if (c == 0) rho[ndxc] = rhoA + 0.2*(r1 - 0.5);
					if (c == 1) rho[ndxc] = rhoB + 0.2*(r1 - 0.5);
					rhosum += rho[ndxc];
					u[ndxc] = 0.0;
					v[ndxc] = 0.0;
				}
				rhoTot[ndx] = rhosum;
				uTot[ndx] = 0.0;
				vTot[ndx] = 0.0;
			}
		}
	}

	//	---------------------------------------
   //	Get the equilibrium distributions:
   //	---------------------------------------

	//exchangeHaloRho();
	fluidForces();
	equilDist();

	//	---------------------------------------
   //	Particle distributions assigned equil.
	// values:
   //	---------------------------------------

   for (int i=0; i<nxyzcn; i++) {
      f[i] = feq[i];
      fprev[i] = f[i];
   }

}



// -------------------------------------------------------------------------
// Take one step forward in the Lattice-Boltzmann simulation:
// {Note: we are using the Shan-Chen multi-component, multi-phase model}
// -------------------------------------------------------------------------

void mcmp::updateMCMP()
{
   //	---------------------------------------
   //	update fprev = f:
   //	---------------------------------------

   for (int i=0; i<nxyzcn; i++) fprev[i] = f[i];

   //	---------------------------------------
   //	execute steps of the LBM:
   //	---------------------------------------

   updateMacros();
   fluidForces();
   equilDist();
   collisionStep();
   streamingStep();
}



// -------------------------------------------------------------------------
// Update macros (density, velocity):
// -------------------------------------------------------------------------

void mcmp::updateMacros()
{
   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {       // loop over grid
            int ndx = i*deli + j*delj + k*delk;
            double sumrho  = 0.0;
            double sumrhou = 0.0;
            double sumrhov = 0.0;
            for (int c=0; c<nc; c++) {      // loop over components
               int ndxc = ndx*nc + c;
               double sum  = 0.0;
               double sumx = 0.0;
               double sumy = 0.0;
               for (int n=0; n<nn; n++) {   // loop over neighbors
                  int ndxcn = ndxc*nn + n;
                  sum  += f[ndxcn];
                  sumx += f[ndxcn]*ex[n];
                  sumy += f[ndxcn]*ey[n];
               }
               rho[ndxc] = sum;
               u[ndxc] = sumx/rho[ndxc];
               v[ndxc] = sumy/rho[ndxc];
               sumrho  += rho[ndxc];
               sumrhou += rho[ndxc]*u[ndxc];
               sumrhov += rho[ndxc]*v[ndxc];
            }
            rhoTot[ndx] = sumrho;
            uTot[ndx] = sumrhou/rhoTot[ndx];
            vTot[ndx] = sumrhov/rhoTot[ndx];
         }
      }
   }
}



// -------------------------------------------------------------------------
// Calculate intermolecular fluid forces:
// -------------------------------------------------------------------------

void mcmp::fluidForces()
{
   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {

            int ndx = i*deli + j*delj + k*delk;
            int ndxA = ndx*nc + 0;
            int ndxB = ndx*nc + 1;

            double GsumxA = 0.0;
            double GsumyA = 0.0;
            double GsumxB = 0.0;
            double GsumyB = 0.0;
            for (int n=0; n<nn; n++) {
               int inbr = nbrIndex(i,int(ex[n]),nx);
         	   int jnbr = nbrIndex(j,int(ey[n]),ny);
               int nbr  = inbr*deli + jnbr*delj + k*delk;
               int nbrA = nbr*nc + 0;
               int nbrB = nbr*nc + 1;
         	   GsumxA += wa[n]*psi(rho[nbrA])*ex[n];
         	   GsumyA += wa[n]*psi(rho[nbrA])*ey[n];
         	   GsumxB += wa[n]*psi(rho[nbrB])*ex[n];
         	   GsumyB += wa[n]*psi(rho[nbrB])*ey[n];
            }

            fx[ndxA] = -G*psi(rho[ndxA])*GsumxB;
            fy[ndxA] = -G*psi(rho[ndxA])*GsumyB;
            fx[ndxB] = -G*psi(rho[ndxB])*GsumxA;
            fy[ndxB] = -G*psi(rho[ndxB])*GsumyA;

         }
      }
   }
}



// -------------------------------------------------------------------------
// Calculate equilibrium distribution, feq:
// -------------------------------------------------------------------------

void mcmp::equilDist()
{
   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {
            int ndx = i*deli + j*delj + k*delk;
            for (int c=0; c<nc; c++) {
               int ndxc = ndx*nc + c;
               // modified velocities due to forces...
               double u_F = uTot[ndx]+ tau*fx[ndxc]/rho[ndxc];
               double v_F = vTot[ndx]+ tau*fy[ndxc]/rho[ndxc];
               double uv2 = u_F*u_F + v_F*v_F;
               for (int n=0; n<nn; n++) {
                  int ndxcn = ndxc*nn + n;
                  double evel = ex[n]*u_F + ey[n]*v_F;
                  feq[ndxcn] = rho[ndxc]*wa[n]*(1.0+3.0*evel+4.5*evel*evel-1.5*uv2);
               }
            }
         }
      }
   }
}



// -------------------------------------------------------------------------
// Collision step:
// -------------------------------------------------------------------------

void mcmp::collisionStep()
{
   for (int i=0; i<nxyzcn; i++) {
      double fdiff = fprev[i] - feq[i];
      omega[i] = -fdiff/tau;
   }
}



// -------------------------------------------------------------------------
// Streaming step:
// -------------------------------------------------------------------------

void mcmp::streamingStep()
{
   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {
            for (int c=0; c<nc; c++) {
               for (int n=0; n<nn; n++) {
                  int ndxcn = (i*deli + j*delj + k*delk)*nc*nn + c*nn + n;
                  int inbr  = nbrIndex(i,int(ex[n]),nx);
                  int jnbr  = nbrIndex(j,int(ey[n]),ny);
                  int nbrcn = (inbr*deli + jnbr*delj + k*delk)*nc*nn + c*nn + n;
                  f[nbrcn] = fprev[ndxcn] + omega[ndxcn];
               }
            }
         }
      }
   }
}



// -------------------------------------------------------------------------
// Potential function of the fluid density:
// -------------------------------------------------------------------------

double mcmp::psi(double rho_in)
{
   return (1.0 - exp(-rho_in));
}



// -------------------------------------------------------------------------
// Alternative potential function of the fluid density:
// -------------------------------------------------------------------------

double mcmp::psi2(double rho_in)
{
	double rho0 = 0.7;
	return rho0*(1.0 - exp(-rho_in/rho0));
}



// -------------------------------------------------------------------------
// neighbor index calculation:
// -------------------------------------------------------------------------

int mcmp::nbrIndex(int i, int di, int ni)
{
   int val = i + di;
   if (val < 1)  val = val + ni;
   if (val > ni) val = val - ni;
   return val;
}



// -------------------------------------------------------------------------
// Print VTK output file:
// -------------------------------------------------------------------------

void mcmp::writeVTKFile(std::string tagname, int tagnum, int iskip, int jskip, int kskip)
{

   // -----------------------------------
   //	Define the file location and name:
   // -----------------------------------

   ofstream outfile;
   std::stringstream filenamecombine;
   filenamecombine << "vtkoutput/" << tagname << "_" << tagnum << ".vtk";
   string filename = filenamecombine.str();
   outfile.open(filename.c_str(), ios::out | ios::app);

   // -----------------------------------
   //	Write the 'vtk' file header:
   // -----------------------------------

   if (rank == 0) {
      string d = "   ";
      outfile << "# vtk DataFile Version 3.1" << endl;
      outfile << "VTK file containing grid data" << endl;
      outfile << "ASCII" << endl;
      outfile << " " << endl;
      outfile << "DATASET STRUCTURED_POINTS" << endl;
      outfile << "DIMENSIONS" << d << NX/iskip << d << NY/jskip << d << NZ/kskip << endl;
      outfile << "ORIGIN " << d << 1 << d << 1 << d << 1 << endl;
      outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0*kskip << endl;
      outfile << " " << endl;
      outfile << "POINT_DATA " << (NX/iskip)*(NY/jskip)*(NZ/kskip) << endl;
      outfile << "SCALARS " << tagname << " float" << endl;
      outfile << "LOOKUP_TABLE default" << endl;
   }

   MPI::COMM_WORLD.Barrier();

   // -----------------------------------
   //	Write the data:
   // NOTE: x-data increases fastest,
   //       then y-data, then z-data
   // -----------------------------------

   int np = MPI::COMM_WORLD.Get_size();    // # of processors

   for (int k=1; k<nz+1; k+=kskip) {
      for (int j=1; j<ny+1; j+=jskip) {
         for (int r=0; r<np; r++) {
            if (r == rank) {
               for (int i=1; i<nx+1; i++) {
                  int ig = i + xOff;
                  if (ig == 0 || ig%iskip == 0) {
                     int ndxA = (i*deli + j*delj + k*delk)*nc + 0;  // print rhoA
                     outfile << fixed << setprecision(3) << rho[ndxA] << endl;
                  }
               }
            }
            MPI::COMM_WORLD.Barrier();
         }
      }
   }

   // -----------------------------------
   //	Close the file:
   // -----------------------------------

   outfile.close();

}
