
# ifndef LBSYSTEM_H
# define LBSYSTEM_H

# include <mpi.h>
# include <string>
# include "../utils/GetPot"

class LBSystem {

    private:

        // basic domain data:

        GetPot input_params;
        int rank,np;
        int nbrL,nbrR;
        int current_step;
        int halo;
        int nn;  // # of nearest neighbors
        int nnn; // # of nearest & next-nearest neighbors
        int nc;
        int nxGlobal;
        int nx,ny,nz;
        int xOff;
        int istr,jstr,kstr;
        int iend,jend,kend;
        double dx,dy,dz;
        double dt;

        // arrays:

        int*** solid;
        int***** nbr;
        double*** rhoTot;
        double*** uTot;
        double*** vTot;
        double**** rho;
        double**** u;
        double**** v;
        double**** fx;
        double**** fy;
        double***** f;
        double***** feq;
        double***** fprev;
        double***** omega;

        // lattice vectors and weights:

        std::string lattice_type;
        double* ex;
        double* ey;
        double* wa;
        double* pa;

        // LB parameters:

        double tau, mu;


    public:

        LBSystem(const GetPot&);
        ~LBSystem();
        void setTimeStep(int step) {current_step = step;}
        void parseInitialCondition();
        void calculateMacros();
        void parseFluidForces();
        void equilibriumDistribution();
        void collisionStep();
        void streamingStep();
        void bounceBack();
        void writeVTKFile(std::string,int,int,int,int);

    private:

        void initializeSolidChannel();
        void initializeSolidFromFile();
        void initializeSingleDroplet();
        void initializeBinaryMixed();
        void initializeBinaryFoam();
        void calculateFluidForcesMCMP();
        void calculateFluidForcesMCMPfoam();
        void addBodyForceXDir();
        void exchangeHaloRho();
        void exchangeHaloDist();
        int nbrIndex(int,int,int);
        double psi(double);
        double psi2(double);

};

# endif  // LBSYSTEM_H
