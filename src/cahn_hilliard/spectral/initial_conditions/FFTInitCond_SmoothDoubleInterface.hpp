

# ifndef FFTSMOOTHDOUBLEINTERFACE_H
# define FFTSMOOTHDOUBLEINTERFACE_H

# include "FFTInitCond_Interface.hpp"

/*
   This is an implementation of the 'FFTInitCond' interface class.
   It initializes a single variable with two smooth, linear interfaces
   defined by a tanh() function.
*/

class FFTSmoothDoubleInterface: public FFTInitCond {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int nx,ny,nz,xOff;
   int deli,delj,delk;
   int sizeR;
   double dx;
   FFTArrays** c;
   vector<double>& eta;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTSmoothDoubleInterface(const Params& p, FFTArrays** cin,
                             vector<double>& etain) :
                             eta(etain)
   {
      c = cin;
      nc = p.nc;
      nx = p.nx;
      ny = p.ny;
      nz = p.nz;
      dx = p.dx;
      deli = (2*(p.nz/2+1))*p.ny;
   	delj = (2*(p.nz/2+1));
   	delk = 1;
      xOff = p.xOff;
      sizeR = p.nx*p.ny*(2*(p.nz/2+1));
      srand(time(NULL)*(p.rank+1));   // set the random seed
   }


   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTSmoothDoubleInterface()
   {
   }

   // -----------------------------
   // Function to calculate i.c.:
   // -----------------------------

   void icFunc()
   {
      for (int n=0; n<nc; n++) {
         for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
               for (int k=0; k<nz; k++) {
                  int ndx = i*deli + j*delj + k*delk;
                  double xi = dx*double(i+xOff);
                  double val = 0.5*(tanh((xi-200)/2.0) + 1.0) + 0.5*(-tanh((xi-100)/2.0) + 1.0);
                  c[n]->setConc(ndx,val);
               }
            }
         }
      }
   }

};

# endif  // FFTSMOOTHDOUBLEINTERFACE_H
