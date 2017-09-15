
# ifndef MCMP_H
# define MCMP_H

# include <mpi.h>
# include <string>
# include <vector>
# include "../utils/GetPot"

class mcmp {

    private:

        // basic domain data:

        GetPot input_params;
        int rank,np;
        int nbrL,nbrR;
        int current_step;
        int halo;
        int nn;  // # of nearest neighbors
        int nc;  // # of components (should be 2)
        int NX,NY,NZ;
        int nx,ny,nz;
        int xOff;
        int istr,jstr,kstr;
        int iend,jend,kend;
        int deli,delj,delk;
        int nxyz,nxyzc,nxyzcn;
        double dx,dy,dz;
        double dt;
        double tau,mu;
        double G;
        std::string lattice_type;
        std::vector<double> rhoTot;
        std::vector<double> uTot;
        std::vector<double> vTot;
        std::vector<double> rho;
        std::vector<double> u,v;
        std::vector<double> fx,fy;
        std::vector<double> f,feq,fprev;
        std::vector<double> omega;
        std::vector<double> ex,ey,wa;

    public:

        mcmp(const GetPot&);
        ~mcmp();
        void setTimeStep(int step) {current_step = step;}
        void initializeMCMP();
        void updateMCMP();
        void writeVTKFile(std::string,int,int,int,int);

    private:

        void updateMacros();
        void fluidForces();
        void equilDist();
        void collisionStep();
        void streamingStep();
        int nbrIndex(int,int,int);
        double psi(double);
        double psi2(double);

};

# endif  // MCMP_H
