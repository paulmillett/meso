
# ifndef PARTICLESFIPI_H
# define PARTICLESFIPI_H

# include <mpi.h>
# include <string>
# include <vector>
# include "Params.hpp"
# include "Sfield.hpp"
# include "Vfield.hpp"

class ParticlesFIPI {

private:

	const Params& p;
	int rank;
	int N;
	int ncell, nncells;
	int ncellx,ncelly,ncellz;
	int current_step;
	double dt, dtover2;
	double dtParticles;
	double rcut, rcut2;
	double pmob;
	double A, scl;
	double cellWidth;
	double cellWidthx,cellWidthy,cellWidthz;
	double Khertz;
	double smallest_rij;
	std::vector <double> r,v,f;
   std::vector <double> mob;
   std::vector <double> rad;
	std::vector <double> box;
	std::vector <double> jam;
	std::vector <int> head, list;
	std::vector <int> cellmap;
	std::vector <int> coord;	

public:

	ParticlesFIPI(const Params&);
	~ParticlesFIPI();
	void initParticles();
   void zeroParticleForces();
   void moveParticles();
	void pairwiseForces();
   void applyBoundaryConditions();
   void writeVTKFile(std::string,int);
	Vfield calcParticleInterfaceForce(const Sfield&, const Vfield&);
	Sfield mapToGrid();
	void calcCapillaryForce(const Sfield&, const Sfield&);
	void setTimeStep(int step) {current_step = step;}

private:

	void setupParticleCells();
	int cellIndex(int,int,int);
	void calcIJInteraction(int,int);

};

# endif  // PARTICLESFIPI_H
