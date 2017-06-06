
# ifndef BCPZONE_H
# define BCPZONE_H

# include "../PFBaseClass.hpp"
# include "../Sfield.hpp"
# include "../Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class BCPZone : public PFBaseClass {

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
   double alpha;
   double wzone;
   double vzone;
   fftw_plan p_forward;
   fftw_plan p_backward;
	fftw_complex* dummy;

public:

   BCPZone(const CommonParams&, const GetPot&);
   ~BCPZone();
   void initPhaseField();
   void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

private:

   void calculateKfields();

};

# endif  // BCPZONE_H
