
# include "ParticlesFIPI_MPI.hpp"
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

ParticlesFIPI_MPI::ParticlesFIPI_MPI(const Params& pin) : p(pin)
{

   //	---------------------------------------
   //	Get input parameters:
   //	---------------------------------------

   Nglobal = p.N;
   NmaxLoc = 2*(Nglobal/p.np);  // maximum number of local particles
   box.reserve(3);
   box[0] = p.Lx;
   box[1] = p.Ly;
   box[2] = p.Lz;
   dt = p.dt;
   pmob = p.pmob;       // particle mobility coefficient
   rcut = p.rcut;       // cut-off radius for particle-particle interactions
   A = p.A;             // magnitude of capillary adhesion
   scl = p.scl;         // scale capillary forces when applied to fluid
   Khertz = p.Khertz;   // Hertz contact force coefficient   
   rcut2 = rcut*rcut;
   dtover2 = dt/2.0;
   borderL = double(p.xOff)*p.dx;
   borderR = borderL + p.nx*p.dx;

   //	---------------------------------------
   //	Determine processor information:
   //	---------------------------------------

   rank = p.rank;

   //	---------------------------------------
   //	Establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<NmaxLoc; i++) {
      rad.push_back(1.0);
      mob.push_back(1.0);
      for (int k=0; k<3; k++) {
         r.push_back(0.0);
         v.push_back(0.0);
         f.push_back(0.0);
         sndbuf0.push_back(0.0);
         sndbuf1.push_back(0.0);
         rcvbuf0.push_back(0.0);
         rcvbuf1.push_back(0.0);
      }
   }

   //	---------------------------------------
   //	Create random seed for r.n.g.:
   //	---------------------------------------

   srand(time(NULL)*(rank+1));

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

ParticlesFIPI_MPI::~ParticlesFIPI_MPI()
{

}



// -------------------------------------------------------------------------
// Initialize particle positions:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::initParticles()
{

   //	---------------------------------------
   // initialize particles in a cubic array:
   //	---------------------------------------

   int nx = int(pow(N,0.334));
   int ny = nx;
   int nz = nx;
   double dr[3] = {box[0]/nx, box[1]/ny, box[2]/nz};
   double ic[3] = {0.0, 0.0, 0.0};
   int parti = 0;
   for (int ix=0; ix<nx; ix++) {
      ic[0]=ix;
      for (int iy=0; iy<ny; iy++) {
         ic[1]=iy;
         for (int iz=0; iz<nz; iz++) {
            ic[2]=iz;
            for (int k=0; k<3; k++) {
               r[parti*3+k] = dr[k]*(double(ic[k])+0.5);
            }
            parti++;
         }
      }
   }
   N = parti - 1;

   //	---------------------------------------
   // initialize particle velocities:
   //	---------------------------------------

   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) {
         v[i*3+k] = 0.0;
      }
   }

   //	---------------------------------------
   // initialize linked-list cells:
   //	---------------------------------------

   setupParticleCells();

}



// -------------------------------------------------------------------------
// Setup the linked-list cell structure.
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::setupParticleCells()
{

   //	---------------------------------------
   // First, how many cells should exist:
   //	---------------------------------------

   cellWidth = rcut;
   ncellx = int(floor(box[0]/cellWidth));
   ncelly = int(floor(box[1]/cellWidth));
   ncellz = int(floor(box[2]/cellWidth));
   if (p.nz == 1) ncellz = 1;
   ncell = ncellx*ncelly*ncellz;

   //	---------------------------------------
   // Second, establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<ncell; i++) {
      head.push_back(-1);
      if (p.nz == 1) nncells = 4;
      if (p.nz  > 1) nncells = 13;
      for (int i=0; i<nncells; i++) {
         cellmap.push_back(0);
      }
   }

   for (int i=0; i<NmaxLoc; i++) {
      list.push_back(0);
   }

   //	---------------------------------------
   // Last, create the cell nabor map:
   //	---------------------------------------

   for (int i=0; i<ncellx; i++) {
      for (int j=0; j<ncelly; j++) {
         for (int k=0; k<ncellz; k++) {

            int imap = cellIndex(i,j,k)*nncells;
            cellmap[imap+0]  = cellIndex(i+1, j , k );
            cellmap[imap+1]  = cellIndex(i+1,j+1, k );
            cellmap[imap+2]  = cellIndex( i ,j+1, k );
            cellmap[imap+3]  = cellIndex(i-1,j+1, k );
            if (nncells == 4) continue; // only do below if 3D cell structure
            cellmap[imap+4]  = cellIndex(i+1, j ,k-1);
            cellmap[imap+5]  = cellIndex(i+1,j+1,k-1);
            cellmap[imap+6]  = cellIndex( i ,j+1,k-1);
            cellmap[imap+7]  = cellIndex(i-1,j+1,k-1);
            cellmap[imap+8]  = cellIndex(i+1, j ,k+1);
            cellmap[imap+9]  = cellIndex(i+1,j+1,k+1);
            cellmap[imap+10] = cellIndex( i ,j+1,k+1);
            cellmap[imap+11] = cellIndex(i-1,j+1,k+1);
            cellmap[imap+12] = cellIndex( i , j ,k+1);

         }
      }
   }

}



// -------------------------------------------------------------------------
// Function to return the linked-list cell index given the
// x- y- z-coordinates of cell:
// -------------------------------------------------------------------------

int ParticlesFIPI_MPI::cellIndex(int i, int j, int k)
{
   if (i < 0) i += ncellx;
   if (i >= ncellx) i -= ncellx;
   if (j < 0) j += ncelly;
   if (j >= ncelly) j -= ncelly;
   if (k < 0) k += ncellz;
   if (k >= ncellz) k -= ncellz;
   return i*ncellz*ncelly + j*ncellz + k;
}



// -------------------------------------------------------------------------
// Re-set particle forces to zero:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::zeroParticleForces()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) f[i*3+k] = 0.0;
   }
}



// -------------------------------------------------------------------------
// Step forward for the particles:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::moveParticles()
{

   //	---------------------------------------
	//	Compute velocities from forces:
	//	---------------------------------------

   for (int i=0; i<N; i++) {
      // for (int k=0; k<3; k++) v[i*3+k] = pmob*f[i*3+k];
      v[i*3+0] = pmob*0.1;
      v[i*3+1] = 0.0;
      v[i*3+2] = 0.0;
   }

   //	---------------------------------------
	//	Update positions:
	//	---------------------------------------

   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] += dt*v[i*3+k];
   }

   //	---------------------------------------
	//	Enforce PBC's:
	//	---------------------------------------

   applyBoundaryConditions();

}



// -------------------------------------------------------------------------
// Pairwise particle forces:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::pairwiseForces()
{

   //	---------------------------------------
   //	First, update linked-list vectors:
   //	---------------------------------------

   for (int i=0; i<ncell; i++) head[i] = -1;

   for (int i=0; i<N; i++) {
      int icell = int(floor(r[i*3+0]/cellWidth))*ncellz*ncelly +
                  int(floor(r[i*3+1]/cellWidth))*ncellz +
                  int(floor(r[i*3+2]/cellWidth));
      list[i] = head[icell];
      head[icell] = i;
   }

   //	---------------------------------------
   //	Then, loop over cells to calculate
   // interactions:
   //	---------------------------------------

   for (int ic=0; ic<ncellx; ic++)
   for (int jc=0; jc<ncelly; jc++)
   for (int kc=0; kc<ncellz; kc++) {

      // loop over particles in current cell:
      int icell = ic*ncellz*ncelly + jc*ncellz + kc;
      int i = head[icell];
      while (i >= 0) {

         // loop over other particles in current cell:
         int j = list[i];
         while (j >= 0) {
            calcIJInteraction(i,j);
            j = list[j];
         }

         // loop over neighboring cells:
         for (int nbor=0; nbor<nncells; nbor++) {
            int jcell = cellmap[icell*nncells + nbor];
            // loop over all particles in jcell:
            int k = head[jcell];
            while (k >= 0) {
               calcIJInteraction(i,k);
               k = list[k];
            }
         }

         // move to next particle in icell:
         i = list[i];

      }

   }

	//	---------------------------------------
   //	Calculate pairwise forces (slow way):
   //	---------------------------------------

   // for (int i=0; i<N; i++) {
   //    for (int j=0; j<i; j++) {
   //       calcIJInteraction(i,j);
   //    }
   // }

}



// -------------------------------------------------------------------------
// Calculate force between particle "i" and "j":
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::calcIJInteraction(int i, int j)
{
   //compute the squared particle distance:
   double dr[3];
   double rij2 = 0.0;
   for (int k=0; k<3; k++) {
      dr[k] = r[i*3+k] - r[j*3+k];
      dr[k] -= round(dr[k]/box[k])*box[k];  // <-- pbc's
      rij2 += dr[k]*dr[k];
   }
   // compute inter-particle forces within cut-off distance:
   if (rij2 <= rcut2) {
      double rij = sqrt(rij2);             // center-to-center dist.
      double s2s = rij - (rad[i]+rad[j]);  // surface-to-surface dist.
      double fij = Khertz*pow((rcut-rij),1.5);    // Hertz contact force
      for (int k=0; k<3; k++) {
         f[i*3+k] += fij*dr[k]/rij;
         f[j*3+k] -= fij*dr[k]/rij;
      }
   }
}



// -------------------------------------------------------------------------
// Enforce boundary conditions (periodic):
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::applyBoundaryConditions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] -= floor(r[i*3+k]/box[k])*box[k];
   }
}



// -------------------------------------------------------------------------
// Calculate the force due to particle-fluid-interface interaction:
// -------------------------------------------------------------------------

Vfield ParticlesFIPI_MPI::calcParticleInterfaceForce(const Sfield& c, const Vfield& gradc)
{

   //	---------------------------------------
   //	Loop over particles, add contributions:
   //	---------------------------------------

   Vfield force_pi(p);  // particle-fluid-interface force Vfield

   for (int i=0; i<N; i++) {
      double cp = c.interpolate(r[i*3+0],r[i*3+1],r[i*3+2]);
      if (abs(cp) <= 0.964) {
         double dis = abs(0.707*log10((1.0+cp)/(1.0-cp)));
         double fip = A*M_PI*0.943*dis;
         double dxx = -cp*gradc.interpolateX(r[i*3+0],r[i*3+1],r[i*3+2]);
         double dyy = -cp*gradc.interpolateY(r[i*3+0],r[i*3+1],r[i*3+2]);
         double dzz = -cp*gradc.interpolateZ(r[i*3+0],r[i*3+1],r[i*3+2]);
         double drr = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
         double nxx = dxx/drr;
         double nyy = dyy/drr;
         double nzz = dzz/drr;
         double fipx = fip*nxx;
         double fipy = fip*nyy;
         double fipz = fip*nzz;
         // add forces to particle force array:
         f[i*3+0] += fipx;
         f[i*3+1] += fipy;
         f[i*3+2] += fipz;
         // add forces to fluid-force array:
         fipx *= -1.0*scl;
         fipy *= -1.0*scl;
         fipz *= -1.0*scl;
         force_pi.addExtrapolation(fipx,fipy,fipz,r[i*3+0],r[i*3+1],r[i*3+2]);
      }
   }

   return force_pi;

}



// -------------------------------------------------------------------------
// Determine which particles are still in current domain, and move
// particles that have left current domain to neighboring domains:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::borderCrossingMPI()
{

   //	---------------------------------------
   //	Loop over particles, check their x-pos.
   //	---------------------------------------

   int Nstay = 0;  // # of particles staying
   int NlveL = 0;  // # of particles leaving left
   int NlveR = 0;  // # of particles leaving right

   for (int i=0; i<N; i++) {
      double xi = r[i*3+0];
      // decide what to do with particle "i":
      if (xi < borderL) {
         // particle has moved to left neighbor:
         sndbuf0[NlveL*3+0] = r[i*3+0];
         sndbuf0[NlveL*3+1] = r[i*3+1];
         sndbuf0[NlveL*3+2] = r[i*3+2];
         NlveL++;
      } else if (xi >= borderR) {
         // particle has moved to right neighbor:
         sndbuf1[NlveR*3+0] = r[i*3+0];
         sndbuf1[NlveR*3+1] = r[i*3+1];
         sndbuf1[NlveR*3+2] = r[i*3+2];
         NlveR++;
      } else {
         // particle is staying in this domain:
         r[Nstay*3+0] = r[i*3+0];
         r[Nstay*3+1] = r[i*3+1];
         r[Nstay*3+2] = r[i*3+2];
         Nstay++;
      }
   }
   N = Nstay;

   //	---------------------------------------
   //	Perform the data communications:
   //	---------------------------------------

   int NentrL = communicateMPI(NlveL,0);
   int NentrR = communicateMPI(NlveR,1);

   //	---------------------------------------
   //	Add immigrants' data to resident data:
   //	---------------------------------------

   for (int i=0; i<NentrL; i++) {
      r[N*3+0] = rcvbuf0[i*3+0];
      r[N*3+1] = rcvbuf0[i*3+1];
      r[N*3+2] = rcvbuf0[i*3+2];
      N++;
   }

   for (int i=0; i<NentrR; i++) {
      r[N*3+0] = rcvbuf1[i*3+0];
      r[N*3+1] = rcvbuf1[i*3+1];
      r[N*3+2] = rcvbuf1[i*3+2];
      N++;
   }

}



// -------------------------------------------------------------------------
// Determine which particles are near the domain borders, and send their
// positions to neighboring domains:
// -------------------------------------------------------------------------

// void ParticlesFIPI_MPI::haloExchangeMPI()
// {
//
//    //	---------------------------------------
//    //	Loop over particles, see if they are
//    // in 'halo' regions:
//    //	---------------------------------------
//
//    int N
//
//    for (int i=0; i<N; i++) {
//
//    }
// }



// -------------------------------------------------------------------------
// Swap array data using MPI commands:
// Input: size = # of particles emigrating in the 'dir' direction
// -------------------------------------------------------------------------

int ParticlesFIPI_MPI::communicateMPI(int size, int dir)
{

}



// -------------------------------------------------------------------------
// Write VTK file for particles:
// -------------------------------------------------------------------------

void ParticlesFIPI_MPI::writeVTKFile(string tagname, int tagnum)
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

   string d = "   ";
   outfile << "# vtk DataFile Version 3.1" << endl;
   outfile << "VTK file containing particle data" << endl;
   outfile << "ASCII" << endl;
   outfile << " " << endl;
   outfile << "DATASET POLYDATA" << endl;
   outfile << " " << endl;
   outfile << "POINTS" << d << N << d << " float" << endl;

   // -----------------------------------
   //	Write the data:
   // NOTE: x-data increases fastest,
   //       then y-data, then z-data
   // -----------------------------------

   for (int i=0; i<N; i++) {
      outfile << fixed << setprecision(3) << r[i*3+0] << d << r[i*3+1] << d << r[i*3+2] << endl;
   }

   // -----------------------------------
   //	Close the file:
   // -----------------------------------

	outfile.close();
}
