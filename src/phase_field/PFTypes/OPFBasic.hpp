
# ifndef OPFBASIC_H
# define OPFBASIC_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class OPFBasic : public PFBaseClass {

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

   OPFBasic(const CommonParams&, const GetPot&);
   ~OPFBasic();
   void initPhaseField();
   void updatePhaseField();
   void outputPhaseField();
   void setTimeStep(int step) {current_step = step;}

private:

   void calculateKfields();

};

# endif  // BCPBASIC_H
