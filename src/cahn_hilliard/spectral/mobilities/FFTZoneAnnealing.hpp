
# ifndef FFTZONEANNEALING_H
# define FFTZONEANNEALING_H

# include "FFTMobilities.hpp"


class FFTZoneAnnealing: public FFTMobilities {

private:

   int nc;
   double vzone;
   double wzone;
   FFTArrays** c;

public:

   FFTZoneAnnealing(const Params&, FFTArrays**);
   ~FFTZoneAnnealing();
   void mobFunc(int,int,int,int,int);

};

# endif  // FFTZONEANNEALING_H
