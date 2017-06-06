

# ifndef NONE_H
# define NONE_H

# include "SourceTerms.hpp"

class None: public SourceTerms {

private:


public:

   None();
   ~None();
   double srcFunc(int,int,int,int,int);

};

# endif  // NONE_H
