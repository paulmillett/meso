
# ifndef TIPSBATHPHIL_H
# define TIPSBATHPHIL_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPSbathPHIL : public PFBaseClass {

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
    double Tbath;
    double Tinit;
    double noiseStr;
    double kappa;
    double nu;
    double gamma;
    double D0;
    double Mweight;
    double Mvolume;
public:

    TIPSbathPHIL(const CommonParams&, const GetPot&);
    ~TIPSbathPHIL();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // TIPSBATHPHIL_H
