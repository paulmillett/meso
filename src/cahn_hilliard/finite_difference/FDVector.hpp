
# ifndef FDVECTOR_H
# define FDVECTOR_H

# include <string>
# include <vector>
# include "Params.hpp"

class FDVector {

private:

	static int instance_count;
	int tag;
	int nxGlobal;
	int xOffset;
	int np;
	int nx;
	int ny;
	int nz;
	int gx;
	int gy;
	int gz;
	int gxyz;
	int deli;
	int delj;
	int delk;
	int nbrL;
	int nbrR;
	int rank;
	double dx,dx2;
	double dy,dy2;
	double dz,dz2;
	std::vector<double> a;

public:

	FDVector(const Params&,double);
	~FDVector();
	double getValue(int) const;
	void setValue(int,double);
	void addValue(int,double);
	void initMeanRandom(double,double);
	void updateBoundaryConditions();
	void mpiBorderExchange();
	void updateExplicitEuler(int,double,double);
	double calculateLaplacian(int) const;
	void writeVTKFile(std::string,int,int,int,int);

};

# endif  // FDVECTOR_H
