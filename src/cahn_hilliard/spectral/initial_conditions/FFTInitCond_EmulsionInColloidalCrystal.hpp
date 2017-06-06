

# ifndef FFTEMULSIONINCOLLOIDALCRYSTAL_H
# define FFTEMULSIONINCOLLOIDALCRYSTAL_H

# include "FFTInitCond_Interface.hpp"

/*
   This is an implementation of the 'FFTInitCond' interface class.
   It initializes a random 2-phase fluid throughout an
   FCC colloidal crystal.
*/

class FFTEmulsionInColloidalCrystal: public FFTInitCond {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int sizeR;
   int nx,ny,nz;
   int NX,xOff;
   int deli,delj,delk;
   double rColl;
   double c0_init_mean;
   double c0_init_noise;
   FFTArrays** c;
   vector<double>& eta;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTEmulsionInColloidalCrystal(const Params& p, FFTArrays** cin,
                                 vector<double>& etain) :
                                 eta(etain)
   {
      c = cin;
      nc = p.nc;
      nx = p.nx;
      ny = p.ny;
      nz = p.nz;
      NX = p.NX;
      deli = (2*(p.nz/2+1))*p.ny;
   	delj = (2*(p.nz/2+1));
   	delk = 1;
      xOff = p.xOff;
      sizeR = p.nx*p.ny*(2*(p.nz/2+1));
      c0_init_mean = p.input_params("CHFFTApp/initial_condition/c0_init_mean",0.0);
      c0_init_noise = p.input_params("CHFFTApp/initial_condition/c0_init_noise",0.0);
      rColl = p.input_params("CHFFTApp/initial_condition/rColl",1.0);

      srand(time(NULL)*p.rank);  // set random seed
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTEmulsionInColloidalCrystal()
   {
   }

   // -----------------------------
   // Function to calculate i.c.:
   // -----------------------------

   void icFunc()
   {

      //	---------------------------------------
      //	define the eta field that describes
      // colloidal crystal positions:
      //	---------------------------------------

      int cx[4]; int cy[4]; int cz[4];
      cx[0] = 0;    cy[0] = 0;    cz[0] = 0;
      cx[1] = NX/2; cy[1] = ny/2; cz[1] = 0;
      cx[2] = NX/2; cy[2] = 0;    cz[2] = nz/2;
      cx[3] = 0;    cy[3] = ny/2; cz[3] = nz/2;

      for (int n=0; n<4; n++) {
         for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
               for (int k=0; k<nz; k++) {
                  int ndx = i*deli + j*delj + k*delk;
                  int iz = k;
                  int iy = j;
                  int ix = i + xOff;
                  int rx = abs(ix - cx[n]);
                  int ry = abs(iy - cy[n]);
                  int rz = abs(iz - cz[n]);
                  if (rx > NX/2) rx = NX - rx;
                  if (ry > ny/2) ry = ny - ry;
                  if (rz > nz/2) rz = nz - rz;
                  double rr = sqrt(rx*rx + ry*ry + rz*rz);
                  if (rr <= rColl) eta[ndx] = 1.0;
               }
            }
         }
      }

      //	---------------------------------------
      //	define the concentration fields
      //	---------------------------------------

      for (int n=0; n<nc; n++) {
         for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
               for (int k=0; k<nz; k++) {
                  int ndx = i*deli + j*delj + k*delk;
                  double r0 = (double)rand()/RAND_MAX;
                  double val = c0_init_mean + c0_init_noise*(r0-0.5);;
                  if (eta[ndx] == 1.0) val = -1.0;
                  c[n]->setConc(ndx,val);
               }
            }
         }
      }

   }

};

# endif  // FFTEMULSIONINCOLLOIDALCRYSTAL_H
