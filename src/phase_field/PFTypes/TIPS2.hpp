
# ifndef TIPS2_H
# define TIPS2_H

# include "../PFBaseClass.hpp"
# include "../Sfield.hpp"
# include "../Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class TIPS2 : public PFBaseClass {

private:

   const CommonParams& p;
	int current_step;
   int nxyz;
	Sfield c;
   Sfield phi;
	Vfield k1;
	Sfield k2;
	Sfield k4;
   double co;
   double M;
   double w;
   double kapc,kapp;
   double phio;
   double prate;
   fftw_plan p_forward;
   fftw_plan p_backward;
	fftw_complex* dummy;

public:

   TIPS2(const CommonParams&, const GetPot&);
   ~TIPS2();
   void initPhaseField();
   void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

private:

   void calculateKfields();

};

# endif  // TIPS2_H
