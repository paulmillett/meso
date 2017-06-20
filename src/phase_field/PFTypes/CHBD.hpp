
# ifndef CHBD_H
# define CHBD_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include "../PFUtils/ParticlesBDCH.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class CHBD: public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nxyz;
    Sfield c1;
    Sfield c2;
    Vfield k1;
    Sfield k2;
    Sfield k4;
    ParticlesBDCH particles;
    double co;
    double M;
    double w;
    double kap;
    fftw_plan p_forward;
    fftw_plan p_backward;
    fftw_complex* dummy;

public:

    CHBD(const CommonParams&, const GetPot&);
    ~CHBD();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

private:

    void calculateKfields();

};

# endif  // CHBD_H
