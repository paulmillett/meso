
# ifndef PARTICLES_H
# define PARTICLES_H

# include <mpi.h>
# include <string>
# include <vector>
# include "initial_conditions/PInitCond_Interface.hpp"
# include "inter_particle_forces/PInterForce_Interface.hpp"
# include "Params.hpp"

class Particles {

private:

	const Params& p;
	int rank;
	int N;
	int current_step;
	double dt;
	double dtover2;
	double Lx;
	double Ly;
	double Lz;
	double Lang;
	double rcut;
	double rcut2;
	double pmob;
	std::vector <double> r,v,f;
	std::vector <double> mass;
	std::vector <double> rad;
	std::vector <double> box;
	PInitCond* icObj;
	PInterForce* fijObj;

public:

	Particles(const Params&);
	~Particles();
	void setTimeStep(int step) {current_step = step;}
	void initParticles();
	void updateParticles();
	void pairwiseForces();
	void updatePositions();
	void velocityHalfKick();
	void applyBoundaryConditions();
	void writeVTKFile(string,int);

};

# endif  // PARTICLES_H
