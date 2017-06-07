
# ifndef PDBASECLASS_H
# define PDBASECLASS_H

# include "../utils/CommonParams.h"
# include "../utils/GetPot"
# include "initial_conditions/PInitCond_Interface.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for particle-dynamics classes in the PD App.
// This class serves as an interface, and contains a factory method.
// It also contains basic elements of particle simulation (e.g. velocity-
// verlet, etc.)
// ---------------------------------------------------------------------

class PDBaseClass {

    public:

        // -------------------------------------------------------------------
        // Factory method that creates objects of sub-classes:
        // -------------------------------------------------------------------

        static PDBaseClass* PDFactory(const CommonParams&, const GetPot&);

        // -------------------------------------------------------------------
        // pure virtual functions:
        // -------------------------------------------------------------------

        virtual void fijFunc(int,int) = 0;

        // -------------------------------------------------------------------
        // common functions:
        // -------------------------------------------------------------------

        PDBaseClass(const CommonParams&, const GetPot&);
        ~PDBaseClass();
        void updateParticles();
        void outputParticles();
        void setTimeStep(int step) {current_step = step;}
        void initParticles();

    protected:

        int N;
        int rank;
        int current_step;
        std::vector <double> r,v,f;
        std::vector <double> rad;
        std::vector <double> box;
        std::vector <double> mass;

    private:

        int ncell, nncells;
        int ncellx,ncelly,ncellz;
        double dt;
        double dtover2;
        double rcut;
        double rcut2;
        double cellWidth;
        double cellWidthx,cellWidthy,cellWidthz;
        std::vector <int> head,list;
        std::vector <int> cellmap;
        PInitCond* icObj;
        void velocityHalfKick();
        void updatePositions();
        void applyBoundaryConditions();
        void pairwiseForces();
        void writeVTKFile(string,int);
        void setupParticleCells();
        int cellIndex(int,int,int);  
  
};

# endif  // PDBASECLASS_H
