
# ifndef CHSYSTEM_H
# define CHSYSTEM_H

# include <mpi.h>
# include <vector>
# include <string>
# include "FDVector.hpp"
# include "chemical_potentials/ChemPot.hpp"
# include "source_terms/SourceTerms.hpp"
# include "mobilities/Mobilities.hpp"
# include "Params.hpp"

class CHSystem {

private:

	const Params& p;
	int current_step;
	int deli,delj,delk;
	double mu_noise;
	string src_type;
 	vector<FDVector> c;
   vector<FDVector> mu;
   ChemPot* muObj;
	SourceTerms* srcObj;
	Mobilities* mobObj;

public:

	CHSystem(const Params&);
	~CHSystem();
   void setTimeStep(int step) {current_step = step;}
	void initializeCahnHilliard();
   void updateCahnHilliard();
	void writeOutputFiles(int);

};

# endif  // CHSYSTEM_H
