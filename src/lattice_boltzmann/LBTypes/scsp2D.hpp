
# ifndef SCSP2D_H
# define SCSP2D_H

# include "../LBBaseClass.hpp"
# include "../LBUtils/LBfluid2D.hpp"
# include "../LBUtils/D2Q9.hpp"


class scsp2D : public LBBaseClass {

private:

	int current_step;
	int deli,delj;
	int nxy,nxyn;
	double tau,mu;
	LBfluid2D fl;
	const D2Q9 ltc;   // struct

	std::vector<double> rho;
	std::vector<double> u,v;
	std::vector<double> fx,fy;
	std::vector<double> f,feq,fprev;
	std::vector<double> omega;
	std::vector<double> ex,ey,wa;

public:

	scsp2D(const CommonParams&, const GetPot&);
	~scsp2D();
	void setTimeStep(int step) {current_step = step;}
	void initLatticeBoltzmann();
	void updateLatticeBoltzmann();
	void outputLatticeBoltzmann();

private:

	double psi(double);
	double psi2(double);

};

# endif  // SCSP2D_H
