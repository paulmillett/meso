

# ifndef FFTHOMOGENEOUSRANDOMNOISE2_H
# define FFTHOMOGENEOUSRANDOMNOISE2_H

# include "FFTInitCond_Interface.hpp"

/*
   This is an implementation of the 'FFTInitCond' interface class.
   It initializes a random 2-phase fluid throughout an
   FCC colloidal crystal.
*/

class FFTHomogeneousRandomNoise2: public FFTInitCond {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int sizeR;
   double c0_init_mean;
   double c0_init_noise;
   double c1_init_mean;
   double c1_init_noise;
   FFTArrays** c;
   vector<double>& eta;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   FFTHomogeneousRandomNoise2(const Params& p, FFTArrays** cin,
                              vector<double>& etain) :
                              eta(etain)
   {
      c = cin;
      nc = p.nc;
      sizeR = p.nx*p.ny*(2*(p.nz/2+1));
      c0_init_mean = p.input_params("CHFFTApp/initial_condition/c0_init_mean",0.3);
      c0_init_noise = p.input_params("CHFFTApp/initial_condition/c0_init_noise",0.1);
      c1_init_mean = p.input_params("CHFFTApp/initial_condition/c1_init_mean",0.3);
      c1_init_noise = p.input_params("CHFFTApp/initial_condition/c1_init_noise",0.1);

      srand(time(NULL)*(p.rank+1));  // set random seed
   }

   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTHomogeneousRandomNoise2()
   {
   }

   // -----------------------------
   // Function to calculate i.c.:
   // -----------------------------

   void icFunc()
   {
      for (int i = 0; i < sizeR; ++i) {
         double r0 = (double)rand()/RAND_MAX;
         double val0 = c0_init_mean + c0_init_noise*(r0-0.5);
         c[0]->setConc(i,val0);
         double r1 = (double)rand()/RAND_MAX;
         double val1 = c1_init_mean + c1_init_noise*(r1-0.5);
         c[1]->setConc(i,val1);
      }
   }

};

# endif  // FFTHOMOGENEOUSRANDOMNOISE2_H
