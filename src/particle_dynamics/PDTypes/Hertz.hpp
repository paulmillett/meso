
# ifndef HERTZ_H
# define HERTZ_H

# include "../PDBaseClass.hpp"

class Hertz : public PDBaseClass {

private:

   double K;
   const CommonParams& p;

public:

   Hertz(const CommonParams&, const GetPot&);
   ~Hertz();
   void initParticles();
   void fijFunc(int,int);

};

# endif  // HERTZ_H
