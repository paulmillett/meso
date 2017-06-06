

# ifndef ZONEANNEALING_H
# define ZONEANNEALING_H

# include "Mobilities.hpp"

class ZoneAnnealing: public Mobilities {

private:

   int nc;
   double vzone;
   double wzone;
   const vector<FDVector>& c;

public:

   ZoneAnnealing(const Params&,const vector<FDVector>&);
   ~ZoneAnnealing();
   double mobFunc(int,int,int,int,int);

};

# endif  // ZONEANNEALING_H
