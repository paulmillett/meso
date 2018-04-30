
# ifndef OPFZone_H
# define OPFZone_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>
# include <cmath>

class OPFZone : public PFBaseClass {

protected:

	const CommonParams& p;
	int current_step;
	int nxyz;
	Sfield c;
	Vfield k1;
	Sfield k2;
	Sfield k4;
	double co;
	double w;
	double wzone;
	double vzone;
	fftw_plan p_forward;
	fftw_plan p_backward;
	fftw_complex* dummy;
	double chiN;
	double chiN2;
	double chiN3;
	double chiN4;
	double c_2;
	double c_4;
	double c_5;
	double c_7;
	double S;
	double singleWellZero;
public:

   OPFZone(const CommonParams&, const GetPot&);
   ~OPFZone();
   void initPhaseField();
   virtual void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

protected:

   void calculateKfields();

};

# endif  // OPFZone_H
