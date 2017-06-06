
# include "ParticlesFIPI.hpp"
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

ParticlesFIPI::ParticlesFIPI(const Params& pin) : p(pin)
{

   //	---------------------------------------
   //	Get input parameters:
   //	---------------------------------------

   N = p.N;
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
   current_step = 0;
   dtParticles = dt/double(p.dtRatio);

   //	---------------------------------------
   //	Determine processor information:
   //	---------------------------------------

   rank = p.rank;

   //	---------------------------------------
   //	Establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<N; i++) {
      rad.push_back(1.0);
      mob.push_back(1.0);
      jam.push_back(0.0);
      coord.push_back(0);
      for (int k=0; k<3; k++) {
         r.push_back(0.0);
         v.push_back(0.0);
         f.push_back(0.0);
      }
   }

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

ParticlesFIPI::~ParticlesFIPI()
{

}



// -------------------------------------------------------------------------
// Initialize particle positions:
// -------------------------------------------------------------------------

void ParticlesFIPI::initParticles()
{

   //	---------------------------------------
   // initialize particles randomly (3D):
   //	---------------------------------------

   // if (p.nz > 1) {
   //    for (int i=0; i<N; i++) {
   //       // assign random position:
   //       r[i*3+0] = (double)rand()/RAND_MAX*box[0];
   //       r[i*3+1] = (double)rand()/RAND_MAX*box[1];
   //       r[i*3+2] = (double)rand()/RAND_MAX*box[2];
   //       // check if it overlaps with existing position:
   //       for (int j=0; j<i; j++) {
   //          double dr[3];
   //          double rij2 = 0.0;
   //          for (int k=0; k<3; k++) {
   //             dr[k] = r[i*3+k] - r[j*3+k];
   //             dr[k] -= round(dr[k]/box[k])*box[k];  // <-- pbc's
   //             rij2 += dr[k]*dr[k];
   //          }
   //          if (rij2 <= rcut2) {
   //             i = i - 1;
   //             continue;  // leave the current loop (j)
   //          }
   //       }
   //    }
   // }

   //	---------------------------------------
   // initialize particles in a cubic array:
   // (3D)
   //	---------------------------------------

   if (p.nz > 1) {
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
   }

   //	---------------------------------------
   // initialize particles in a square array:
   // (2D)
   //	---------------------------------------

   if (p.nz == 1) {
      int nx = int(pow(N,0.5));
      int ny = nx;
      int nz = 1;
      double drx = box[0]/nx;
      double dry = box[1]/ny;
      int parti = 0;
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            r[parti*3+0] = drx*(ix+0.5);
            r[parti*3+1] = dry*(iy+0.5);
            r[parti*3+2] = 0.0;
            parti++;
         }
      }
   }

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

void ParticlesFIPI::setupParticleCells()
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
   cellWidthx = box[0]/double(ncellx);
   cellWidthy = box[1]/double(ncelly);
   cellWidthz = box[2]/double(ncellz);

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

   for (int i=0; i<N; i++) {
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

int ParticlesFIPI::cellIndex(int i, int j, int k)
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

void ParticlesFIPI::zeroParticleForces()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) f[i*3+k] = 0.0;
   }
}



// -------------------------------------------------------------------------
// Step forward for the particles:
// -------------------------------------------------------------------------

void ParticlesFIPI::moveParticles()
{

   //	---------------------------------------
	//	Compute particle mobilities:
	//	---------------------------------------

   for (int i=0; i<N; i++) {
      mob[i] = pmob * (0.5*(1 - tanh(3*(jam[i]-1.2)))); 
      //mob[i] = pmob;
   }

   //	---------------------------------------
	//	Compute velocities from forces:
	//	---------------------------------------

   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) v[i*3+k] = mob[i]*f[i*3+k];
   }

   //	---------------------------------------
	//	Update positions:
	//	---------------------------------------

   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] += dtParticles*v[i*3+k];
   }

   //	---------------------------------------
	//	Enforce PBC's:
	//	---------------------------------------

   applyBoundaryConditions();

}



// -------------------------------------------------------------------------
// Pairwise particle forces:
// -------------------------------------------------------------------------

void ParticlesFIPI::pairwiseForces()
{

   //	---------------------------------------
   //	zero out the coordinations and 'jam':
   //	---------------------------------------

   for (int i=0; i<N; i++) {
      coord[i] = 0;
      jam[i] = 0.0;
   }

   //	---------------------------------------
   //	First, update linked-list vectors:
   //	---------------------------------------

   for (int i=0; i<ncell; i++) head[i] = -1;

   for (int i=0; i<N; i++) {
      int icell = int(floor(r[i*3+0]/cellWidthx))*ncellz*ncelly +
                  int(floor(r[i*3+1]/cellWidthy))*ncellz +
                  int(floor(r[i*3+2]/cellWidthz));
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

void ParticlesFIPI::calcIJInteraction(int i, int j)
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
      //double fij = 2/pow(rij,5);   // polynomial 'soft' force
      for (int k=0; k<3; k++) {
         f[i*3+k] += fij*dr[k]/rij;
         f[j*3+k] -= fij*dr[k]/rij;
      }
      coord[i]++;
      coord[j]++;
      jam[i] += rcut - rij;
      jam[j] += rcut - rij;
   }
}



// -------------------------------------------------------------------------
// Enforce boundary conditions (periodic):
// -------------------------------------------------------------------------

void ParticlesFIPI::applyBoundaryConditions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] -= floor(r[i*3+k]/box[k])*box[k];
   }
}



// -------------------------------------------------------------------------
// Calculate the force due to particle-fluid-interface interaction:
// -------------------------------------------------------------------------

Vfield ParticlesFIPI::calcParticleInterfaceForce(const Sfield& c, const Vfield& gradc)
{

   //	---------------------------------------
   //	Loop over particles, add contributions:
   //	---------------------------------------

   Vfield force_pi(p);  // particle-fluid-interface force Vfield

   double cbound = 0.964;
   for (int i=0; i<N; i++) {
      double cp = c.interpolate(r[i*3+0],r[i*3+1],r[i*3+2]);
      if (abs(cp) <= cbound) {

         // define particle-interface force:
         //double dis = abs(0.707*log10((1.0+cp)/(1.0-cp))); // Botto's equation
         //double fip = A*M_PI*0.943*dis;                    // Botto's equation
         double fip = A*M_PI*0.943*(abs(cp)*(cbound-abs(cp)));  // my equation
         //double fip = A*M_PI*0.943*(abs(cp));                 // my equation #2

         // calculate direction to interface center:
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
// Map particle positions onto an Sfield grid:
// -------------------------------------------------------------------------

Sfield ParticlesFIPI::mapToGrid()
{
   // width of region surrounding particles to probe:
   int wdth = 7;
   // loop over particles:
   Sfield eta(p);
   for (int pp=0; pp<N; pp++) {
      // particle's position:
      double x = r[pp*3+0];
      double y = r[pp*3+1];
      double z = r[pp*3+2];
      // grid point nearest particle (rounded down):
      int x0 = int(floor(x)/p.dx);
      int y0 = int(floor(y)/p.dy);
      int z0 = int(floor(z)/p.dz);
      // loop over region near particle:
      for (int i=0; i<wdth; i++) {
         int ii = x0 - (wdth/2 - 1) + i;
         if (ii < 0) ii += p.NX;
         if (ii > p.NX-1) ii -= p.NX;
         for (int j=0; j<wdth; j++) {
            int jj = y0 - (wdth/2 - 1) + j;
            if (jj < 0) jj += p.NY;
            if (jj > p.NY-1) jj -= p.NY;
            for (int k=0; k<wdth; k++) {
               int kk = z0 - (wdth/2 - 1) + k;
               if (kk < 0) kk += p.NZ;
               if (kk > p.NZ-1) kk -= p.NZ;
               if (p.NZ == 1) kk = 0;
               // calculate distance to point:
               double rx = x - double(ii);
               double ry = y - double(jj);
               double rz = z - double(kk);
               rx -= round(rx/(p.NX*p.dx))*p.NX*p.dx;
               ry -= round(ry/(p.NY*p.dy))*p.NY*p.dy;
               rz -= round(rz/(p.NZ*p.dz))*p.NZ*p.dz;
               double r2 = rx*rx + ry*ry + rz*rz;
               // assign spread function to grid:
               double val = 1.0 - r2/6.25; //1.0/(1.0 + 4.0*r2);
               if (val < 0.0) val = 0.0;
               double val0 = creal(eta.getValue(ii*p.nz*p.ny + jj*p.nz + kk));
               eta.setValue(ii*p.nz*p.ny + jj*p.nz + kk, max(val,val0));
            }
         }
      }
   }
   return eta;
}



// -------------------------------------------------------------------------
// Calculate capillary forces on particles due to fluid-fluid interface:
// -------------------------------------------------------------------------

void ParticlesFIPI::calcCapillaryForce(const Sfield& eta, const Sfield& c)
{
   // width of region surrounding particles to probe:
   int wdth = 7;
   // loop over particles:
   for (int pp=0; pp<N; pp++) {
      // particle's position:
      double x = r[pp*3+0];
      double y = r[pp*3+1];
      double z = r[pp*3+2];
      // grid point nearest particle (rounded down):
      int x0 = int(floor(x)/p.dx);
      int y0 = int(floor(y)/p.dy);
      int z0 = int(floor(z)/p.dz);
      // loop over region near particle:
      for (int i=0; i<wdth; i++) {
         int ii = x0 - (wdth/2 - 1) + i;
         if (ii < 0) ii += p.NX;
         if (ii > p.NX-1) ii -= p.NX;
         for (int j=0; j<wdth; j++) {
            int jj = y0 - (wdth/2 - 1) + j;
            if (jj < 0) jj += p.NY;
            if (jj > p.NY-1) jj -= p.NY;
            for (int k=0; k<wdth; k++) {
               int kk = z0 - (wdth/2 - 1) + k;
               if (kk < 0) kk += p.NZ;
               if (kk > p.NZ-1) kk -= p.NZ;
               if (p.NZ == 1) kk = 0;
               double et1 = creal(eta.getValue(ii*p.nz*p.ny + jj*p.nz + kk));
               if (et1 >= 0.0) {
                  // calculate distance to point:
                  double rx = double(ii) - x;
                  double ry = double(jj) - y;
                  double rz = double(kk) - z;
                  rx -= round(rx/(p.NX*p.dx))*p.NX*p.dx;
                  ry -= round(ry/(p.NY*p.dy))*p.NY*p.dy;
                  rz -= round(rz/(p.NZ*p.dz))*p.NZ*p.dz;
                  double r2 = rx*rx + ry*ry + rz*rz;
                  // assign spread function to grid:
                  double et2 = 1.0 - r2/6.25;
                  if (et2 < 0.0) et2 = 0.0;
                  double cijk = creal(c.getValue(ii*p.nz*p.ny + jj*p.nz + kk));
                  double cint = pow(1-cijk*cijk,8);
                  double fint = et2*cint;
                  f[pp*3+0] += A*fint*rx;
                  f[pp*3+1] += A*fint*ry;
                  f[pp*3+2] += A*fint*rz;
               }
            }
         }
      }
   }
}



// -------------------------------------------------------------------------
// Write VTK file for particles:
// -------------------------------------------------------------------------

void ParticlesFIPI::writeVTKFile(string tagname, int tagnum)
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

   outfile << " " << endl;
   outfile << "POINT_DATA " << N << endl;
   outfile << "SCALARS Coord float" << endl;
   outfile << "LOOKUP_TABLE deault" << endl;

   for (int i=0; i<N; i++) {
      //outfile << fixed << setprecision(1) << double(coord[i]) << endl;
      outfile << fixed << setprecision(3) << jam[i] << endl;
   }

   // -----------------------------------
   //	Close the file:
   // -----------------------------------

	outfile.close();
}
