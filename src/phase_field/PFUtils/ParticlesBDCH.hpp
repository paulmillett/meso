
# ifndef PARTICLESBDCH_H
# define PARTICLESBDCH_H

# include "../../particle_dynamics/PDParticles.hpp"
# include "Sfield.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is a particle class that inherits from PDParticles and adds
// specific features for coupling particle data to phase-field data.  In
// particlar, features for calculating capillary forces between CH
// interfaces and particles.
// ---------------------------------------------------------------------

class ParticlesBDCH : public PDParticles {

    public:

        ParticlesBDCH(const CommonParams&, const GetPot&);
        ~ParticlesBDCH();
        Sfield mapToGrid();
        void calcCapillaryForce(const Sfield&, const Sfield&, const Sfield&);

    protected:

        string chbdType;
        bool thinFilm;
        const CommonParams& p;
        double cap_str;
        double eps,n;
        int thickness;
        std::vector <double> fcap;
        std::vector <double> fcapSum;
        void auxiliaryForces();
        bool isMyParticle(double,int);
        bool inMyDomain(int);

};

# endif  // PARTICLESBDCH_H
