
# ifndef TIPS3_H
# define TIPS3_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include <complex.h>
# include <fftw3-mpi.h>
# include <mpi.h>


class TIPS3 : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nxyz;
    Sfield c;
    Vfield k1;
    Sfield k2;
    Sfield k4;
    double co;
    double M;
    double N;
    double kT;
    double chi;
    double kap;
    double A,B;
    fftw_plan p_forward;
    fftw_plan p_backward;
    fftw_complex* dummy;

public:

    TIPS3(const CommonParams&, const GetPot&);
    ~TIPS3();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

private:

    void calculateKfields();

};

# endif  // TIPS3_H
