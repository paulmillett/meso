# ifndef OPFZoneTempFD_H
# define OPFZoneTempFD_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class OPFZoneTempFD : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nx,ny,nz;
    int deli,delj,delk;
    int nxyz;
    SfieldFD c;
    double co;
    double M;
    double noiseStr;
	double wzone;
	double vzone;
	double S;
	double Tmax;
	double Tmin;
	double N;
	double alpha;
	double beta;
	int templating;
	int templateSpacing;

public:

    OPFZoneTempFD(const CommonParams&, const GetPot&);
    ~OPFZoneTempFD();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // OPFZoneTempFD_H
