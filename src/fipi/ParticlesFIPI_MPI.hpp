
# ifndef PARTICLESFIPI_MPI_H
# define PARTICLESFIPI_MPI_H

# include <mpi.h>
# include <string>
# include <vector>
# include "Params.hpp"
# include "Sfield.hpp"
# include "Vfield.hpp"

class ParticlesFIPI_MPI {

private:

	const Params& p;
	int rank;
	int N,Nglobal,NmaxLoc;
	int ncell, nncells;
	int ncellx,ncelly,ncellz;
	int current_step;
	double borderL,borderR;
	double dt, dtover2;
	double rcut, rcut2;
	double pmob;
	double A, scl;
	double cellWidth;
	double Khertz;
	std::vector <double> r,v,f;
   std::vector <double> mob;
   std::vector <double> rad;
	std::vector <double> box;
	std::vector <double> sndbuf0,sndbuf1;
	std::vector <double> rcvbuf0,rcvbuf1;
	std::vector <int> head, list;
	std::vector <int> cellmap;

public:

	ParticlesFIPI_MPI(const Params&);
	~ParticlesFIPI_MPI();
	void initParticles();
   void zeroParticleForces();
   void moveParticles();
	void pairwiseForces();
   void applyBoundaryConditions();
   void writeVTKFile(std::string,int);
	Vfield calcParticleInterfaceForce(const Sfield&, const Vfield&);

private:

	void setupParticleCells();
	int cellIndex(int,int,int);
	void calcIJInteraction(int,int);
	void borderCrossingMPI();
	int communicateMPI(int,int);

};

# endif  // PARTICLESFIPI_MPI_H
