

# ifndef RANDOM_H
# define RANDOM_H

# include "PInitCond_Interface.hpp"

/*
   This is an implementation of the 'PInitCond' interface class.
   It initializes a set of particles at random positions in the domain.
   It requires particles to be separated by a minimum distance; consequently,
   if you try to use a number of particles that results in a particle volume
   fraction above about 10% it will take FOREVER to initialize. So only use
   this initial condition for particle volume fractions <= 10%.
   */

class Random : public PInitCond {

    // -------------------------------------------------------------------------
    // Private class members:
    // -------------------------------------------------------------------------

    private:

        int N;
        double Lx, Ly, Lz;
        double vscl;
        double pradii;
        vector<double>& r;
        vector<double>& v;
        vector<double>& rad;

        // -------------------------------------------------------------------------
        // Public class methods:
        // -------------------------------------------------------------------------

    public:

        // -----------------------------
        // Constructor...
        // -----------------------------

        Random(const GetPot& p, vector<double>& rin,
                vector<double>& vin, vector<double>& radin) :
            r(rin), v(vin), rad(radin)
    {
        N = p("BDApp/N",1);
        Lx = p("BDApp/Lx",5.0);
        Ly = p("BDApp/Ly",5.0);
        Lz = p("BDApp/Lz",5.0);
        vscl = p("BDApp/initial_condition/vscl",0.0);
        pradii = p("BDApp/initial_condition/pradii",1.0);
    }


        // -----------------------------
        // Destructor...
        // -----------------------------

        ~Random()
        {
        }

        // -----------------------------
        // Function to calculate i.c.:
        // -----------------------------

        void icFunc()
        {
            // initialize particle positions at random locations in the domain
            for (int i=0; i<N; i++)
            {
                bool tooClose = true;
                double r1,r2,r3;
                while (tooClose)
                {
                    tooClose = false;

                    // get a random position
                    r1 = (double)rand()/RAND_MAX*Lx;
                    r2 = (double)rand()/RAND_MAX*Ly;
                    r3 = (double)rand()/RAND_MAX*Lz;

                    // assign position
                    r[i*3 + 0] = r1;
                    r[i*3 + 1] = r2;
                    r[i*3 + 2] = r3;

                    // check to see if random position is too close to other positions
                    for (int k = 0; k<i; k++)
                    {
                        double drx = calc_separation_pbc(r[3*i + 0],r[3*k + 0],Lx);
                        double dry = calc_separation_pbc(r[3*i + 1],r[3*k + 1],Ly);
                        double drz = calc_separation_pbc(r[3*i + 2],r[3*k + 2],Lz);
                        double rij = sqrt(drx*drx+dry*dry+drz*drz);
                        if (rij < 3.0*pradii) // assumes all particles have same radius
                        {
                            tooClose = true;
                            break;
                        }
                    }
                }
            }

            // initialize particle radii:
            for (int i=0; i<N; i++) {
                rad[i] = pradii;
            }

            // initialize particle velocities:
            double vsum[3] = {0.0, 0.0, 0.0};
            for (int i=0; i<N; i++) {
                for (int k=0; k<3; k++) {
                    double r = (double)rand()/RAND_MAX;
                    v[i*3+k] = vscl*(r - 0.5);
                    vsum[k] += v[i*3+k];
                }
            }

            // zero the total momentum:
            for (int i=0; i<N; i++) {
                for (int k=0; k<3; k++) {
                    v[i*3+k] -= vsum[k]/N;
                }
            }

        }

        double calc_separation_pbc(double ri, double rj, double L)
        {
            double rij = abs(ri - rj);
            if (rij <= L/2.0)
            {
                return rij;
            }
            else
                return L - rij;
        }

};

# endif  // RANDOM_H
