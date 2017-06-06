
# ifndef FFTSIMPLEDOUBLEWELL_H
# define FFTSIMPLEDOUBLEWELL_H

# include "FFTChemPot.hpp"


class FFTSimpleDoubleWell: public FFTChemPot {

private:

   int nc;
   double w;
   double kap;
   FFTArrays** c;

public:

   FFTSimpleDoubleWell(const Params&, FFTArrays**);
   ~FFTSimpleDoubleWell();
   void muFunc(int,int,int,int,int);

};

# endif  // FFTSIMPLEDOUBLEWELL_H
