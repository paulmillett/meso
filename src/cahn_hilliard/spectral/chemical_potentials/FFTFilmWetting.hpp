
# ifndef FFTFILMWETTING_H
# define FFTFILMWETTING_H

# include "FFTChemPot.hpp"


class FFTFilmWetting: public FFTChemPot {

private:

   int nc;
   int nz;
   double w;
   double kap;
   double bias_str;
   FFTArrays** c;

public:

   FFTFilmWetting(const Params&, FFTArrays**);
   ~FFTFilmWetting();
   void muFunc(int,int,int,int,int);

};

# endif  // FFTFILMWETTING_H
