
# ifndef FIPISYSTEM_H
# define FIPISYSTEM_H

# include <mpi.h>
# include "Sfield.hpp"
# include "Vfield.hpp"
# include "ParticlesFIPI.hpp"
# include "Params.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class FIPISystem {

private:

	const Params& p;
	int current_step;
   int nxyz;
	double visc;
	Sfield c;
	Vfield u;
	Vfield k1;
	Sfield k2;
	Sfield k4;
	ParticlesFIPI particles;
   fftw_plan p_forward;
   fftw_plan p_backward;
	fftw_complex* dummy;

public:

	FIPISystem(const Params&);
	~FIPISystem();
   void setTimeStep(int step) {current_step = step;}
	void initializeFIPI();
   void updateFIPI();
	void writeOutputFiles(int);

};

# endif  // FIPISYSTEM_H
