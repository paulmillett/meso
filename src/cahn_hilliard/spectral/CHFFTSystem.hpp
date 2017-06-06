
# ifndef CHFFTSYSTEM_H
# define CHFFTSYSTEM_H

# include <mpi.h>
# include <vector>
# include <string>
# include "FFTArrays.hpp"
# include "initial_conditions/FFTInitCond_Interface.hpp"
# include "chemical_potentials/FFTChemPot_Interface.hpp"
# include "source_terms/FFTSourceTerms_Interface.hpp"
# include "mobilities/FFTMobilities_Interface.hpp"
# include "Params.hpp"

class CHFFTSystem {

private:

	const Params& p;
	int current_step;
	int deli,delj,delk;
   int sizeR;
	double  c_noise;
	double mu_noise;
	std::vector<double> eta;
	string src_type;
	string mob_type;
	bool* mob_flag;
	FFTArrays** c;
	FFTInitCond* icObj;
 	FFTChemPot* muObj;
	FFTSourceTerms* srcObj;
	FFTMobilities* mobObj;

public:

	CHFFTSystem(const Params&);
	~CHFFTSystem();
   void setTimeStep(int step) {current_step = step;}
	void initializeCahnHilliard();
   void updateCahnHilliard();
	void writeOutputFiles(int);

};

# endif  // CHFFTSYSTEM_H
