
# include "Sfield.hpp"
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

Sfield::Sfield(const Params& pin) : p(pin)
{

   //	---------------------------------------
   //	Unpack some of the 'params' data:
   //	---------------------------------------

   rank = p.rank;
   NX = p.NX;
   NY = p.NY;
   NZ = p.NZ;
   nx = p.nx;
   ny = p.ny;
   nz = p.nz;
   nxyz = nx*ny*nz;
   xOffset = p.xOff;

   //	---------------------------------------
   //	Create array and initialize:
   //	---------------------------------------

   ptrdiff_t locsize, locnx, offx;
   // I can change the below line...!
   locsize = fftw_mpi_local_size_3d(p.NX,p.NY,p.NZ,MPI_COMM_WORLD,&locnx,&offx);
   a = fftw_alloc_complex(locsize);
   for (int i=0; i<nxyz; i++) a[i] = 0.0;
   in_fspace = false;

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

Sfield::~Sfield()
{
   fftw_free(a);
}



// -------------------------------------------------------------------------
// Setter, Getter, & Adder:
// -------------------------------------------------------------------------

void Sfield::setValue(int i, double val)
{
   a[i] = val;
}

void Sfield::setValue(int i, fftw_complex val)
{
   a[i] = val;
}

void Sfield::addValue(int i, double val)
{
   a[i] += val;
}

fftw_complex Sfield::getValue(int i) const
{
   return a[i];
}



// -------------------------------------------------------------------------
// Reset field to zero:
// -------------------------------------------------------------------------

void Sfield::resetSfield(std::string space)
{
   for (int i=0; i<nxyz; i++) a[i] = 0.0+0.0i;
   if (space == "f") in_fspace = true;
   if (space == "r") in_fspace = false;
}



// -------------------------------------------------------------------------
// Interpolate field value at a point between grid nodes:
// {note: see wikipedia page on trilinear interpolation}
// -------------------------------------------------------------------------

double Sfield::interpolate(double x, double y, double z) const
{
   int x0 = int(floor(x)/p.dx);
   int x1 = x0 + 1;
   if (x1 >= nx) x1 = 0;
   int y0 = int(floor(y)/p.dy);
   int y1 = y0 + 1;
   if (y1 >= ny) y1 = 0;
   int z0 = int(floor(z)/p.dz);
   int z1 = z0 + 1;
   if (z1 >= nz) z1 = 0;
   double xd = x - double(x0);
   double yd = y - double(y0);
   double zd = z - double(z0);

   double val = creal(a[x0*nz*ny + y0*nz + z0])*(1.0-xd)*(1.0-yd)*(1.0-zd);
   val += creal(a[x1*nz*ny + y0*nz + z0])*xd*(1.0-yd)*(1.0-zd);
   val += creal(a[x0*nz*ny + y1*nz + z0])*(1.0-xd)*yd*(1.0-zd);
   val += creal(a[x0*nz*ny + y0*nz + z1])*(1.0-xd)*(1.0-yd)*zd;
   val += creal(a[x1*nz*ny + y0*nz + z1])*xd*(1.0-yd)*zd;
   val += creal(a[x0*nz*ny + y1*nz + z1])*(1.0-xd)*yd*zd;
   val += creal(a[x1*nz*ny + y1*nz + z0])*xd*yd*(1.0-zd);
   val += creal(a[x1*nz*ny + y1*nz + z1])*xd*yd*zd;
   return val;

}



// -------------------------------------------------------------------------
// Extrapolate point to grid evenly:
// -------------------------------------------------------------------------

void Sfield::extrapolatePointToGrid(double x, double y, double z, int wdth)
{
   int x0 = int(floor(x)/p.dx);
   int y0 = int(floor(y)/p.dy);
   int z0 = int(floor(z)/p.dz);

   if (wdth == 1) {
      int x1 = x0 + 1;
      if (x1 >= nx) x1 = 0;
      int y1 = y0 + 1;
      if (y1 >= ny) y1 = 0;
      int z1 = z0 + 1;
      if (z1 >= nz) z1 = 0;
      a[x0*nz*ny + y0*nz + z0] = 1.0;
      a[x1*nz*ny + y0*nz + z0] = 1.0;
      a[x0*nz*ny + y1*nz + z0] = 1.0;
      a[x0*nz*ny + y0*nz + z1] = 1.0;
      a[x1*nz*ny + y1*nz + z0] = 1.0;
      a[x0*nz*ny + y1*nz + z1] = 1.0;
      a[x1*nz*ny + y0*nz + z1] = 1.0;
      a[x1*nz*ny + y1*nz + z1] = 1.0;
   }

   if (wdth > 1) {
      for (int i=0; i<wdth; i++) {
         int ii = x0 - (wdth/2 - 1) + i;
         if (ii < 0) ii += NX;
         if (ii > NX-1) ii -= NX;
         for (int j=0; j<wdth; j++) {
            int jj = y0 - (wdth/2 - 1) + j;
            if (jj < 0) jj += NY;
            if (jj > NY-1) jj -= NY;
            for (int k=0; k<wdth; k++) {
               int kk = z0 - (wdth/2 - 1) + k;
               if (kk < 0) kk += NZ;
               if (kk > NZ-1) kk -= NZ;
               if (NZ == 1) kk = 0;
               // calculate distance to point:
               double rx = x - double(ii);
               double ry = y - double(jj);
               double rz = z - double(kk);
               rx -= round(rx/(NX*p.dx))*NX*p.dx;
               ry -= round(ry/(NY*p.dy))*NY*p.dy;
               rz -= round(rz/(NZ*p.dz))*NZ*p.dz;
               double r2 = rx*rx + ry*ry; // + rz*rz;
               // assign spread function to grid:
               double val = 1.0 - r2/16.0; //1.0/(1.0 + 4.0*r2);
               double val0 = 0.0; //creal(a[ii*nz*ny + jj*nz + kk]);
               a[ii*nz*ny + jj*nz + kk] = val; //max(val,val0);
            }
         }
      }
   }

}



// -------------------------------------------------------------------------
// FFTw transforms:
// -------------------------------------------------------------------------

void Sfield::fft(const fftw_plan& p_forward)
{
   fftw_mpi_execute_dft(p_forward,a,a);
   in_fspace = true;
}

void Sfield::ifft(const fftw_plan& p_backward)
{
   fftw_mpi_execute_dft(p_backward,a,a);
   for (int i=0; i<nxyz; i++) a[i] /= NX*NY*NZ;
   in_fspace = false;
}



// -------------------------------------------------------------------------
// Write array values to 'vtk' file:
// -------------------------------------------------------------------------

void Sfield::writeVTKFile(std::string tagname, int tagnum,
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
   	outfile << "DIMENSIONS" << d << NX/iskip << d << NY/jskip << d << NZ/kskip << endl;
   	outfile << "ORIGIN " << d << 0 << d << 0 << d << 0 << endl;
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

   for (int k=0; k<nz; k+=kskip) {
      for (int j=0; j<ny; j+=jskip) {
         for (int r=0; r<np; r++) {
            if (r == rank) {
               for (int i=0; i<nx; i++) {
                  int ig = i + xOffset;
                  if (ig == 0 || ig%iskip == 0) {
                     outfile << fixed << setprecision(3)
                                      << creal(a[i*ny*nz + j*nz + k]) << endl;
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



// -------------------------------------------------------------------------
// Compound Assignment Operators:
// -------------------------------------------------------------------------

Sfield& Sfield::operator+=(const Sfield& rhs)
{
   for (int i=0; i<nxyz; i++) a[i] += rhs.a[i];
   return *this;
}

Sfield& Sfield::operator+=(double val)
{
   for (int i=0; i<nxyz; i++) a[i] += val;
   return *this;
}

Sfield& Sfield::operator-=(const Sfield& rhs)
{
   for (int i=0; i<nxyz; i++) a[i] -= rhs.a[i];
   return *this;
}

Sfield& Sfield::operator-=(double val)
{
   for (int i=0; i<nxyz; i++) a[i] -= val;
   return *this;
}

Sfield& Sfield::operator*=(const Sfield& rhs)
{
   for (int i=0; i<nxyz; i++) a[i] *= rhs.a[i];
   return *this;
}

Sfield& Sfield::operator*=(double val)
{
   for (int i=0; i<nxyz; i++) a[i] *= val;
   return *this;
}

Sfield& Sfield::operator*=(fftw_complex val)
{
   for (int i=0; i<nxyz; i++) a[i] *= val;
   return *this;
}

Sfield& Sfield::operator/=(const Sfield& rhs)
{
   for (int i=0; i<nxyz; i++) {
      if (rhs.a[i] != 0.0) a[i] /= rhs.a[i];
      else a[i] = 0.0;
   }
   return *this;
}

Sfield& Sfield::operator/=(double val)
{
   for (int i=0; i<nxyz; i++) {
      if (val != 0.0) a[i] /= val;
      else a[i] = 0.0;
   }
   return *this;
}

Sfield& Sfield::operator=(const Sfield& rhs)
{
   for (int i=0; i<nxyz; i++) a[i] = rhs.a[i];
   in_fspace = rhs.in_fspace;
   return *this;
}



// -------------------------------------------------------------------------
// Binary Operators:
// -------------------------------------------------------------------------

Sfield Sfield::operator+(const Sfield& rhs) const
{
   Sfield result(p);
   result = *this;
   result += rhs;
   return result;
}

Sfield Sfield::operator+(double val) const
{
   Sfield result(p);
   result = *this;
   result += val;
   return result;
}

Sfield Sfield::operator-(const Sfield& rhs) const
{
   Sfield result(p);
   result = *this;
   result -= rhs;
   return result;
}

Sfield Sfield::operator-(double val) const
{
   Sfield result(p);
   result = *this;
   result -= val;
   return result;
}

Sfield Sfield::operator*(const Sfield& rhs) const
{
   Sfield result(p);
   result = *this;
   result *= rhs;
   return result;
}

Sfield Sfield::operator*(double val) const
{
   Sfield result(p);
   result = *this;
   result *= val;
   return result;
}

Sfield Sfield::operator*(fftw_complex val) const
{
   Sfield result(p);
   result = *this;
   result *= val;
   return result;
}

Sfield Sfield::operator/(const Sfield& rhs) const
{
   Sfield result(p);
   result = *this;
   result /= rhs;
   return result;
}

Sfield Sfield::operator/(double val) const
{
   Sfield result(p);
   result = *this;
   result /= val;
   return result;
}



// -------------------------------------------------------------------------
// Some non-member methods:
// -------------------------------------------------------------------------

const Sfield operator+(double b, const Sfield& a)
{
   return a+b;
}

const Sfield operator-(double b, const Sfield& a)
{
   return a*(-1.0) + b;
}

const Sfield operator*(double b, const Sfield& a)
{
   return a*b;
}
