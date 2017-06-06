

# ifndef FFTDIBLOCKCOPOLYMER_H
# define FFTDIBLOCKCOPOLYMER_H

# include "FFTSourceTerms.hpp"

class FFTDiblockCopolymer: public FFTSourceTerms {

private:

   int nc;
   double dt;
   double A0;
   double c0_init_mean;
   FFTArrays** c;

public:

   FFTDiblockCopolymer(const Params&, FFTArrays**);
   ~FFTDiblockCopolymer();
   void srcFunc(int,int);

};

# endif  // FFTDIBLOCKCOPOLYMER_H
