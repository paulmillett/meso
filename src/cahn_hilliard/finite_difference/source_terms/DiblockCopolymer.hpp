

# ifndef DIBLOCKCOPOLYMER_H
# define DIBLOCKCOPOLYMER_H

# include "SourceTerms.hpp"

class DiblockCopolymer: public SourceTerms {

private:

   int nc;
   double A0;
   double c0_init_mean;
   const vector<FDVector>& c;

public:

   DiblockCopolymer(const Params&,const vector<FDVector>&);
   ~DiblockCopolymer();
   double srcFunc(int,int,int,int,int);

};

# endif  // DIBLOCKCOPOLYMER_H
