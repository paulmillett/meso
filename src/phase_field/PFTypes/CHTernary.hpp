
# ifndef CHTERNARY_H
# define CHTERNARY_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class CHTernary : public PFBaseClass {

private:

   const CommonParams& p;
	int current_step;
   int nxyz;
	Sfield c1;
   Sfield c2;
   Sfield c3;
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

   CHTernary(const CommonParams&, const GetPot&);
   ~CHTernary();
   void initPhaseField();
   void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

private:

   void calculateKfields();

};

# endif  // CHTERNARY_H
