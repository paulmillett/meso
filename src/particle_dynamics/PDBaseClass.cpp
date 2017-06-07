
# include "PDBaseClass.hpp"
# include <string>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
# include <iostream>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PDTypes/Hertz.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PDBaseClass* PDBaseClass::PDFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

   // -----------------------------------
   // identify the requested object:
   // -----------------------------------

   string pd_type = input_params("PDApp/type","Hertz");

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   if (pd_type == "Hertz") return new Hertz(p,input_params);

}



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

PDBaseClass::PDBaseClass(const CommonParams& p, const GetPot& input_params)
{

   cout << "Hello from PDBaseClass Ctor" << endl;

   //	---------------------------------------
   //	Get parameters from 'CommonParams':
   //	---------------------------------------

   box.reserve(3);
   box[0] = p.LX;
   box[1] = p.LY;
   box[2] = p.LZ;
   dt = p.dt;
   rank = p.rank;
   dtover2 = dt/2.0;
   current_step = 0;
   if (p.nz == 1) flag2D = true;
   if (p.nz  > 1) flag2D = false;

   //	---------------------------------------
   //	Get other parameters from 'GetPot':
   //	---------------------------------------

   N = input_params("PDApp/N",1);              // # of particles
   rcut = input_params("PDApp/rcut",4.0);      // cut-off radius for particle-particle interactions
   rcut2 = rcut*rcut;

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
   // initialize linked-list cells:
   //	---------------------------------------

   setupParticleCells();

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

PDBaseClass::~PDBaseClass()
{

}



// -------------------------------------------------------------------------
// Updater:
// -------------------------------------------------------------------------

void PDBaseClass::updateParticles()
{
   velocityHalfKick();
   updatePositions();
   applyBoundaryConditions();
	pairwiseForces();
   velocityHalfKick();
}



// -------------------------------------------------------------------------
// Update the particle positions by the time step:
// -------------------------------------------------------------------------

void PDBaseClass::updatePositions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] += dt*v[i*3+k];
   }
}



// -------------------------------------------------------------------------
// Update particle velocities by half the time step:
// -------------------------------------------------------------------------

void PDBaseClass::velocityHalfKick()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) v[i*3+k] += dtover2*f[i*3+k]/mass[i];
   }
}



// -------------------------------------------------------------------------
// Enforce boundary conditions (periodic):
// -------------------------------------------------------------------------

void PDBaseClass::applyBoundaryConditions()
{
   for (int i=0; i<N; i++) {
      for (int k=0; k<3; k++) r[i*3+k] -= floor(r[i*3+k]/box[k])*box[k];
   }
}



// -------------------------------------------------------------------------
// Pairwise particle forces:
// -------------------------------------------------------------------------

void PDBaseClass::pairwiseForces()
{

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
            fijFunc(i,j);
            j = list[j];
         }

         // loop over neighboring cells:
         for (int nbor=0; nbor<nncells; nbor++) {
            int jcell = cellmap[icell*nncells + nbor];
            // loop over all particles in jcell:
            int k = head[jcell];
            while (k >= 0) {
               fijFunc(i,k);
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
// Outputer:
// -------------------------------------------------------------------------

void PDBaseClass::outputParticles()
{
   writeVTKFile("particles",current_step);
}



// -------------------------------------------------------------------------
// Write VTK file for particles:
// -------------------------------------------------------------------------

void PDBaseClass::writeVTKFile(string tagname, int tagnum)
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



// -------------------------------------------------------------------------
// Setup the linked-list cell structure.
// -------------------------------------------------------------------------

void PDBaseClass::setupParticleCells()
{

   //	---------------------------------------
   // First, how many cells should exist:
   //	---------------------------------------

   cellWidth = rcut;
   ncellx = int(floor(box[0]/cellWidth));
   ncelly = int(floor(box[1]/cellWidth));
   ncellz = int(floor(box[2]/cellWidth));
   if (flag2D) ncellz = 1;
   ncell = ncellx*ncelly*ncellz;
   cellWidthx = box[0]/double(ncellx);
   cellWidthy = box[1]/double(ncelly);
   cellWidthz = box[2]/double(ncellz);

   //	---------------------------------------
   // Second, establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<ncell; i++) {
      head.push_back(-1);
      if (flag2D)  nncells = 4;
      if (!flag2D) nncells = 13;
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

int PDBaseClass::cellIndex(int i, int j, int k)
{
   if (i < 0) i += ncellx;
   if (i >= ncellx) i -= ncellx;
   if (j < 0) j += ncelly;
   if (j >= ncelly) j -= ncelly;
   if (k < 0) k += ncellz;
   if (k >= ncellz) k -= ncellz;
   return i*ncellz*ncelly + j*ncellz + k;
}
