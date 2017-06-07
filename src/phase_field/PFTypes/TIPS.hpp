
# ifndef TIPS_H
# define TIPS_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class TIPS: public PFBaseClass {

private:

   const CommonParams& p;
	int current_step;
   int nxyz;
	Sfield c;
	Vfield k1;
	Sfield k2;
	Sfield k4;
   double co;
   double M;
   double w;
   double kap;
   fftw_plan p_forward;
   fftw_plan p_backward;
	fftw_complex* dummy;

public:

   TIPS(const CommonParams&, const GetPot&);
   ~TIPS();
   void initPhaseField();
   void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

private:

   void calculateKfields();

};

# endif  // TIPS_H
