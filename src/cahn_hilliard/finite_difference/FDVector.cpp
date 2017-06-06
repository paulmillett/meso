
# include "FDVector.hpp"
# include <mpi.h>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// static variable initialization:
// -------------------------------------------------------------------------

int FDVector::instance_count = 0;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FDVector::FDVector(const Params& p, double val)
{

   //	---------------------------------------
   //	Unpack some of the 'params' data:
   //	---------------------------------------

   np = p.np;
   rank = p.rank;
   nbrL = p.nbrL;
   nbrR = p.nbrR;
   nxGlobal = p.NX;
   nx = p.nx;
   ny = p.ny;
   nz = p.nz;
   dx = p.dx;
   dy = p.dy;
   dz = p.dz;
   xOffset = p.xOff;
   dx2 = dx*dx;
   dy2 = dy*dy;
   dz2 = dz*dz;

   //	---------------------------------------
   //	Establish array dimensions:
   //	---------------------------------------

   instance_count++;
   tag = instance_count;   // array identifier
   gx = nx + 2;            // local x-dim. + ghost nodes
   gy = ny + 2;            // local y-dim. + ghost nodes
   gz = nz + 2;            // local z-dim. + ghost nodes
   gxyz = gx*gy*gz;        // total vector size
   deli = gz*gy;           // index offset for neighbors in x-dim.
   delj = gz;              // index offset for neighbors in y-dim.
   delk = 1;               // index offset for neighbors in z-dim.

   //	---------------------------------------
   //	Create vector:
   //	---------------------------------------

   for (int i=0; i<gxyz; i++) {
      a.push_back(val);
   }

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FDVector::~FDVector()
{
}


// -------------------------------------------------------------------------
// Getter:
// -------------------------------------------------------------------------

double FDVector::getValue(int i) const
{
   return a[i];
}



// -------------------------------------------------------------------------
// Setter:
// -------------------------------------------------------------------------

void FDVector::setValue(int i, double val)
{
   a[i] = val;
}



// -------------------------------------------------------------------------
// Adder:
// -------------------------------------------------------------------------

void FDVector::addValue(int i, double val)
{
   a[i] += val;
}



// -------------------------------------------------------------------------
// Initialize array with random values about 'mean' scaled by 'noise':
// -------------------------------------------------------------------------

void FDVector::initMeanRandom(double mean, double noise)
{
   srand(time(NULL)*rank);
   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         for (int k=1; k<nz+1; k++) {
            int ndx = k*delk + j*delj + i*deli;
            double r = (double)rand()/RAND_MAX;
   			a[ndx] = mean + noise*(r-0.5);
         }
		}
	}
}



// -------------------------------------------------------------------------
// Calculate "Laplacian" using given coordinates:
// -------------------------------------------------------------------------

double FDVector::calculateLaplacian(int i) const
{
   double lapx = (a[i+deli] + a[i-deli] - 2*a[i])/(dx2);
   double lapy = (a[i+delj] + a[i-delj] - 2*a[i])/(dy2);
   double lapz = (a[i+delk] + a[i-delk] - 2*a[i])/(dz2);
	return lapx + lapy + lapz;
}



// -------------------------------------------------------------------------
// Update array value using explicit Euler method:
// -------------------------------------------------------------------------

void FDVector::updateExplicitEuler(int i, double dt, double dadt)
{
	a[i] += dt*dadt;
}



// -------------------------------------------------------------------------
// Periodic boundary conditions:
// -------------------------------------------------------------------------

void FDVector::updateBoundaryConditions()
{

   // -----------------------------------
   // Set boundary conditions (x-dir.)
   // -----------------------------------

   mpiBorderExchange();

   // -----------------------------------
   // Set boundary conditions (y-dir.)
   // -----------------------------------

   for (int i=1; i<nx+1; i++) {
      for (int k=1; k<nz+1; k++) {
         a[k + 0*delj + i*deli] = a[k + ny*delj + i*deli];
         a[k + (ny+1)*delj + i*deli] = a[k + 1*delj + i*deli];
      }
	}

   // -----------------------------------
   // Set boundary conditions (z-dir.)
   // -----------------------------------

   for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
         a[0 + j*delj + i*deli] = a[nz + j*delj + i*deli];
         a[(nz+1) + j*delj + i*deli] = a[1 + j*delj + i*deli];
      }
	}

}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors:
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void FDVector::mpiBorderExchange()
{

   MPI::Status status;
   int size = deli;          // size of y-z face that will be sent

   // -----------------------------------
   // Send to left, Recv from right
   // -----------------------------------

   int stamp = tag*10;       // stamp is a unique int for each comm.
   int ondx = 1*deli;        // out index
   int indx = (nx+1)*deli;   // in  index

   MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrL,stamp,
                            &a[indx],size,MPI::DOUBLE,nbrR,stamp,status);

   // -----------------------------------
   // Send to right, Recv from left
   // -----------------------------------

   stamp += 1;               // update the stamp
   ondx = nx*deli;           // out index
   indx = 0;                 // in  index

   MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrR,stamp,
                            &a[indx],size,MPI::DOUBLE,nbrL,stamp,status);

}



// -------------------------------------------------------------------------
// Write array values to 'vtk' file:
// -------------------------------------------------------------------------

void FDVector::writeVTKFile(string tagname, int tagnum,
                            int iskip, int jskip, int kskip)
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
   	outfile << "DIMENSIONS" << d << nxGlobal/iskip << d << ny/jskip << d << nz/kskip << endl;
   	outfile << "ORIGIN " << d << 1 << d << 1 << d << 1 << endl;
   	outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0*kskip << endl;
   	outfile << " " << endl;
   	outfile << "POINT_DATA " << (nxGlobal/iskip)*(ny/jskip)*(nz/kskip) << endl;
   	outfile << "SCALARS " << tagname << " float" << endl;
   	outfile << "LOOKUP_TABLE default" << endl;
   }

   MPI::COMM_WORLD.Barrier();

   // -----------------------------------
   //	Write the data (in VTK format):
   // NOTE: x-data increases fastest,
   //       then y-data, then z-data
   // -----------------------------------

   for (int k=1; k<nz+1; k+=kskip) {
      for (int j=1; j<ny+1; j+=jskip) {
         for (int r=0; r<np; r++) {
            if (r == rank) {
               for (int i=1; i<nx+1; i++) {
                  int ig = i + xOffset;
                  if (ig == 0 || ig%iskip == 0) {
                     int ndx = k*delk + j*delj + i*deli;
                     outfile << fixed << setprecision(3) << a[ndx] << endl;
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
