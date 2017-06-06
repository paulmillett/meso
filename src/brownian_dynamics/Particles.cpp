
# include "Particles.hpp"
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

Particles::Particles(const Params& pin) : p(pin)
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
   rcut = p.input_params("BDApp/inter_particle_forces/rcut",4.0);
   rcut2 = rcut*rcut;
   dtover2 = dt/2.0;

   //	---------------------------------------
   //	Determine processor information:
   //	---------------------------------------

   rank = p.rank;

   //	---------------------------------------
   //	Establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<N; i++) {
      rad.push_back(1.0);
      mass.push_back(1.0);
      for (int k=0; k<3; k++) {
         r.push_back(0.0);
         v.push_back(0.0);
         f.push_back(0.0);
      }
   }

   //	---------------------------------------
   //	Create random seed for r.n.g.:
   //	---------------------------------------

   srand(time(NULL)*(rank+1));

   //	---------------------------------------
   //	Create objects for init. cond. &
   // inter-particle forces:
   //	---------------------------------------

   icObj = PInitCond::PInitCondFactory(p,r,v,rad);
   fijObj = PInterForce::PInterForceFactory(p);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

Particles::~Particles()
{

}



// -------------------------------------------------------------------------
// Initialize particle positions:
// -------------------------------------------------------------------------

void Particles::initParticles()
{
   icObj->icFunc();
}



// -------------------------------------------------------------------------
// Step forward for the particles (velocity-verlet algorithm):
// -------------------------------------------------------------------------

void Particles::updateParticles()
{
   velocityHalfKick();
   updatePositions();
   applyBoundaryConditions();
	pairwiseForces();
   velocityHalfKick();
}



// -------------------------------------------------------------------------
// Pairwise particle forces:
// -------------------------------------------------------------------------

void Particles::pairwiseForces()
{

	//	---------------------------------------
	//	Zero all forces:
	//	---------------------------------------

   double dr[3];

   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) f[i*3+k] = 0.0;
   }

	//	---------------------------------------
   //	Calculate pairwise forces:
   //	---------------------------------------

   for (int i=0; i<N; i++) {
      for (int j=0; j<i; j++) {
         // compute the squared particle distance:
         double rij2 = 0.0;
         for (int k=0; k<3; k++) {
            dr[k] = r[i*3+k] - r[j*3+k];
            dr[k] -= round(dr[k]/box[k])*box[k];   // <-- pbc's
            rij2 += dr[k]*dr[k];
         }
         // compute inter-particle forces within cut-off distance:
         if (rij2 <= rcut2) {
            double rij = sqrt(rij2);             // center-to-center dist.
            double s2s = rij - (rad[i]+rad[j]);  // surface-to-surface dist.
            double fij = fijObj->fijFunc(rij,s2s);
            for (int k=0; k<3; k++) {
               f[i*3+k] += fij*dr[k]/rij;
               f[j*3+k] -= fij*dr[k]/rij;
            }
         }
      }
   }

}



// -------------------------------------------------------------------------
// Update the particle positions by the time step:
// -------------------------------------------------------------------------

void Particles::updatePositions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] += dt*v[i*3+k];
   }
}



// -------------------------------------------------------------------------
// Update particle velocities by half the time step:
// -------------------------------------------------------------------------

void Particles::velocityHalfKick()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) v[i*3+k] += dtover2*f[i*3+k]/mass[i];
   }
}



// -------------------------------------------------------------------------
// Enforce boundary conditions (periodic):
// -------------------------------------------------------------------------

void Particles::applyBoundaryConditions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] -= floor(r[i*3+k]/box[k])*box[k];
   }
}



// -------------------------------------------------------------------------
// Write VTK file for particles:
// -------------------------------------------------------------------------

void Particles::writeVTKFile(string tagname, int tagnum)
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
