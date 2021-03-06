
# ifndef CHBD_H
# define CHBD_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include "../PFUtils/ParticlesBDCH.hpp"
# include "../../utils/rand.hpp"
# include <complex.h>
# include <fftw3-mpi.h>
# include <string>

class CHBD: public PFBaseClass {

    protected:

        const CommonParams& p;
        int current_step;
        int nxyz;
        int part_step_skip;
        int numberOfParticles; //for turning off output when N=0
        Sfield c1;
        Sfield c2;
        Sfield cp;
        Vfield k1;
        Sfield k2;
        Sfield k4;
        Sfield * kz; // for applying E-field in z-dir
        Rand rng; // Mersenne Twister random number generator
        ParticlesBDCH particles;
        double co;
        double M;
        double w;
        double kap;
        double eCH;
        double noiseStr;
        string pfType;
        fftw_plan p_forward;
        fftw_plan p_backward;
        fftw_complex* dummy;

    public:

        CHBD(const CommonParams&, const GetPot&);
        ~CHBD();
        virtual void initPhaseField();
        virtual void updatePhaseField();
        virtual void outputPhaseField();
        virtual void setTimeStep(int step) {current_step = step;}

    protected:

        void calculateKfields();

};

# endif  // CHBD_H
