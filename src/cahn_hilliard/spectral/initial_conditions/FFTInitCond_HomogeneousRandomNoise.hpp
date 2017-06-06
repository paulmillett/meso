

# ifndef FFTHOMOGENEOUSRANDOMNOISE_H
# define FFTHOMOGENEOUSRANDOMNOISE_H

# include "FFTInitCond_Interface.hpp"

/*
   This is an implementation of the 'FFTInitCond' interface class.
   It initializes a single variable with random values averaged with
   c0_init_mean and variation c0_init_noise.
*/

class FFTHomogeneousRandomNoise: public FFTInitCond {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int nc;
   int sizeR;
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

   FFTHomogeneousRandomNoise(const Params& p, FFTArrays** cin,
                             vector<double>& etain) :
                             eta(etain)
   {
      c = cin;
      nc = p.nc;
      sizeR = p.nx*p.ny*(2*(p.nz/2+1));
      c0_init_mean = p.input_params("CHFFTApp/initial_condition/c0_init_mean",0.5);
      c0_init_noise = p.input_params("CHFFTApp/initial_condition/c0_init_noise",0.2);

      srand(time(NULL)*(p.rank+1));   // set the random seed
   }


   // -----------------------------
   // Destructor...
   // -----------------------------

   ~FFTHomogeneousRandomNoise()
   {
   }

   // -----------------------------
   // Function to calculate i.c.:
   // -----------------------------

   void icFunc()
   {
      for (int i = 0; i < sizeR; ++i) {
         double r = (double)rand()/RAND_MAX;
         double val = c0_init_mean + c0_init_noise*(r-0.5);
         c[0]->setConc(i,val);
      }
   }

};

# endif  // FFTHOMOGENEOUSRANDOMNOISE_H
