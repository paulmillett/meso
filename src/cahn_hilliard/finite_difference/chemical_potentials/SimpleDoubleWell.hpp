
# ifndef SIMPLEDOUBLEWELL_H
# define SIMPLEDOUBLEWELL_H

# include "ChemPot.hpp"

class SimpleDoubleWell: public ChemPot {

private:

   int nc;
   double w;
   double kap;
   const vector<FDVector>& c;

public:

   SimpleDoubleWell(const Params&,const vector<FDVector>&);
   ~SimpleDoubleWell();
   double muFunc(int,int,int,int,int);

};

# endif  // SIMPLEDOUBLEWELL_H
