
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

    private:

        const CommonParams& p;
        double cap_str;
        std::vector <double> fcap;
        void auxiliaryForces();

};

# endif  // PARTICLESBDCH_H
