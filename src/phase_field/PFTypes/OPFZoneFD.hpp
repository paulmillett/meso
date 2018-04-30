
# ifndef OPFZoneFD_H
# define OPFZoneFD_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class OPFZoneFD : public PFBaseClass {

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

    OPFZoneFD(const CommonParams&, const GetPot&);
    ~OPFZoneFD();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // OPFZoneFD_H
