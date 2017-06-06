
# include "FFTArrays.hpp"
# include <mpi.h>
# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FFTArrays::FFTArrays(const Params& p)
{

   //	---------------------------------------
   //	Unpack some of the 'params' data:
   //	---------------------------------------

   rank = p.rank;
   nxGlobal = p.NX;
	ny = p.NY;
   nz = p.NZ;
	dx = p.dx;
	dy = p.dy;
   dz = p.dz;
   dt = p.dt;
   int n0 = p.NX;
   int n1 = p.NY;
   int n2 = p.NZ;

   //	---------------------------------------
   //	Initialize FFTw and set array size:
   //	---------------------------------------

   fftw_mpi_init();
   locsize = fftw_mpi_local_size_3d(n0,n1,n2/2+1,MPI_COMM_WORLD,&locnx,&offx);
   c = fftw_alloc_real(2*locsize);
	mu = fftw_alloc_real(2*locsize);
   mob = fftw_alloc_real(2*locsize);
   gradmu = fftw_alloc_real(2*locsize);
   Mgradmu = fftw_alloc_real(2*locsize);
   cFFT = fftw_alloc_complex(locsize);
	muFFT = fftw_alloc_complex(locsize);
   gradmuFFT = fftw_alloc_complex(locsize);
   MgradmuFFT = fftw_alloc_complex(locsize);
   sizeR = 2*locsize;   // size of real arrays
   sizeC = locsize;     // size of compl arrays
   nx = locnx;          // local x-dim.
   xOffset = offx;      // offset for x-dim.
   i1 = I;              // imaginary number

   //	---------------------------------------
   //	Create plans for FFTw transforms:
   //	---------------------------------------

   plan1 = fftw_mpi_plan_dft_r2c_3d(n0,n1,n2,c,cFFT,MPI_COMM_WORLD,FFTW_MEASURE);
	plan2 = fftw_mpi_plan_dft_r2c_3d(n0,n1,n2,mu,muFFT,MPI_COMM_WORLD,FFTW_MEASURE);
   plan3 = fftw_mpi_plan_dft_c2r_3d(n0,n1,n2,cFFT,c,MPI_COMM_WORLD,FFTW_MEASURE);
   plan4 = fftw_mpi_plan_dft_c2r_3d(n0,n1,n2,gradmuFFT,gradmu,MPI_COMM_WORLD,FFTW_MEASURE);
   plan5 = fftw_mpi_plan_dft_r2c_3d(n0,n1,n2,Mgradmu,MgradmuFFT,MPI_COMM_WORLD,FFTW_MEASURE);

   //	---------------------------------------
   //	Create the k^2 and k^4 arrays:
   //	---------------------------------------

   k1 = fftw_alloc_real(locsize);
   k2 = fftw_alloc_real(locsize);
	k4 = fftw_alloc_real(locsize);

   int nx2 = n0/2;
	int ny2 = n1/2;
	int nz2 = n2/2;
	double dkx = (2*M_PI)/(n0*1.0);
	double dky = (2*M_PI)/(n1*1.0);
	double dkz = (2*M_PI)/(n2*1.0);

	for (int i=0; i<locnx; i++) {
		for (int j=0; j<n1; j++) {
			for (int k=0; k<n2/2+1; k++) {
 				int ii = i + offx;
 				int jj = j;
 				int kk = k;
 				double kx = ii*dkx;
 				double ky = jj*dky;
 				double kz = kk*dkz;
 				if (ii > nx2) kx = (ii-n0)*dkx;
 				if (jj > ny2) ky = (jj-n1)*dky;
 				double kx2 = kx*kx;
 				double ky2 = ky*ky;
 				double kz2 = kz*kz;
  				k2[i*n1*(n2/2+1)+j*(n2/2+1)+k] = kx2 + ky2 + kz2;
  				k4[i*n1*(n2/2+1)+j*(n2/2+1)+k] = (kx2 + ky2 + kz2)*(kx2 + ky2 + kz2);
            k1[i*n1*(n2/2+1)+j*(n2/2+1)+k] = sqrt(kx2+ky2+kz2);
 			}
		}
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FFTArrays::~FFTArrays()
{
   fftw_destroy_plan(plan1);
   fftw_destroy_plan(plan2);
	fftw_destroy_plan(plan3);
   fftw_destroy_plan(plan4);
   fftw_destroy_plan(plan5);
}



// -------------------------------------------------------------------------
// Getter for concentration:
// -------------------------------------------------------------------------

double FFTArrays::getConc(int i) const
{
	return c[i];
}



// -------------------------------------------------------------------------
// Setter for concentration:
// -------------------------------------------------------------------------

void FFTArrays::setConc(int i, double val)
{
	c[i] = val;
}



// -------------------------------------------------------------------------
// Add to concentration:
// -------------------------------------------------------------------------

void FFTArrays::addtoConc(int i, double val)
{
	c[i] += val;
}



// -------------------------------------------------------------------------
// Setter for chemical potential:
// -------------------------------------------------------------------------

void FFTArrays::setMu(int i, double val)
{
   mu[i] = val;
}



// -------------------------------------------------------------------------
// Add to chemical potential:
// -------------------------------------------------------------------------

void FFTArrays::addtoMu(int i, double val)
{
	mu[i] += val;
}



// -------------------------------------------------------------------------
// Setter for chemical mobility:
// -------------------------------------------------------------------------

void FFTArrays::setMob(int i, double val)
{
   mob[i] = val;
}



// -------------------------------------------------------------------------
// Getter for concentration:
// -------------------------------------------------------------------------

double FFTArrays::getMob(int i) const
{
	return mob[i];
}



// -------------------------------------------------------------------------
// Getter for real array size:
// -------------------------------------------------------------------------

int FFTArrays::getNxLocal()
{
	return nx;
}



// -------------------------------------------------------------------------
// Getter for real array size:
// -------------------------------------------------------------------------

int FFTArrays::getRealSize()
{
	return sizeR;
}



// -------------------------------------------------------------------------
// Getter for x-offset:
// -------------------------------------------------------------------------

int FFTArrays::getXoffset()
{
	return xOffset;
}



// -------------------------------------------------------------------------
// Perform forward transform:
// -------------------------------------------------------------------------

void FFTArrays::forwardFFT()
{
   fftw_execute(plan1);
	fftw_execute(plan2);
}



// -------------------------------------------------------------------------
// Perform backward transform:
// -------------------------------------------------------------------------

void FFTArrays::backwardFFT()
{
   fftw_execute(plan3);
   for (int i=0; i<sizeR; i++) c[i] /= nxGlobal*ny*nz;
}



// -------------------------------------------------------------------------
// Update 'c' in Fourier space:
// {Note: we follow the approach specified by LQ Chen's spectral CH model
//  given in: Physical Review E (1999) vol. 60, p. 3564}
// -------------------------------------------------------------------------

void FFTArrays::updateConcFFT(bool isomobility)
{

   // -----------------------------------
   //	If mobility is ISOTROPIC:
   // -----------------------------------

   if (isomobility == true) {
      for (int i=0; i<sizeC; i++) {
         cFFT[i] = (cFFT[i] - dt*k2[i]*muFFT[i])/(1.0 + dt*k4[i]);
      }
   }

   // -----------------------------------
   //	If mobility is ANISOTROPIC:
   // {Note: this takes a few more steps}
   // -----------------------------------

   if (isomobility == false) {

      // --------------------------------
      //	calculate gradmuFFT:
      // --------------------------------

      for (int i=0; i<sizeC; i++) {
         gradmuFFT[i] = i1*k1[i]*(muFFT[i] + k2[i]*cFFT[i]);
      }

      // --------------------------------
      //	backward transform 'gradmuFFT':
      // --------------------------------

      fftw_execute(plan4);
      for (int i=0; i<sizeR; i++) gradmu[i] /= nxGlobal*ny*nz;

      // --------------------------------
      //	calculate Mgradmu, then MgradmuFFT:
      // --------------------------------

      for (int i=0; i<sizeR; i++) Mgradmu[i] = mob[i]*gradmu[i];
      fftw_execute(plan5);

      // --------------------------------
      //	update conc. in Fourier space:
      // --------------------------------

      double A = 0.8;
      for (int i=0; i<sizeC; i++) {
         cFFT[i] = cFFT[i] + (dt*i1*k1[i]*MgradmuFFT[i])/(1.0 + A*dt*k4[i]);
      }

   }

}



// -------------------------------------------------------------------------
// Write array values to 'vtk' file:
// -------------------------------------------------------------------------

void FFTArrays::writeVTKFile(string tagname, int tagnum,
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
   //	Write the data:
   // NOTE: x-data increases fastest,
   //       then y-data, then z-data
   // -----------------------------------

   int np = MPI::COMM_WORLD.Get_size();    // # of processors

   for (int k=0; k<nz; k+=kskip) {
      for (int j=0; j<ny; j+=jskip) {
         for (int r=0; r<np; r++) {
            if (r == rank) {
               for (int i=0; i<nx; i++) {
                  int ig = i + xOffset;
                  if (ig == 0 || ig%iskip == 0) {
                     outfile << fixed << setprecision(3)
                                      << mu[(i*ny + j) * (2*(nz/2+1)) + k] << endl;
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
