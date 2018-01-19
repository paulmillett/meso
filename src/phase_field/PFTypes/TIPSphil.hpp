
# ifndef TIPSPHIL_H
# define TIPSPHIL_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPSphil : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nx,ny,nz;
    int deli,delj,delk;
    int nxyz;
    int outAnalysisInterval;
    int numAnalysisOutputs;
    SfieldFD c;
    double co;
    double M;
    double N;
    double alpha;
    double beta;
    double kap;
    double A;
    double Tinit;
    double Tbath;
    double noiseStr;
    double nu;
    double gamma;
    double D0;
    double cHeat;
    double Mweight;
    double Mvolume;
public:

    TIPSphil(const CommonParams&, const GetPot&);
    ~TIPSphil();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // TIPSPHIL_H
