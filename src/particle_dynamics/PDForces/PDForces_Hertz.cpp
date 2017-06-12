
# include "PDForces_Hertz.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

Hertz::Hertz(const GetPot& p, vector<double>& rin, vector<double>& vin,
             vector<double>& fin, vector<double>& radin) :
             r(rin), v(vin), f(fin), rad(radin)
{
   K = p("PDApp/inter_particle_forces/K",0.1);
   box[0] = p("PDApp/Lx",5.0);
   box[1] = p("PDApp/Ly",5.0);
   box[2] = p("PDApp/Lz",5.0);
   rcut = p("PDApp/inter_particle_forces/rcut",4.0);
   rcut2 = rcut*rcut;
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

Hertz::~Hertz()
{
}



// -------------------------------------------------------------------------
// Function to calculate fij:
// -------------------------------------------------------------------------

void Hertz::fijFunc(int i, int j)
{
    //compute the squared particle distance:
    double dr[3];
    double rij2 = 0.0;
    for (int k=0; k<3; k++) {
       dr[k] = r[i*3+k] - r[j*3+k];
       dr[k] -= round(dr[k]/box[k])*box[k];  // <-- pbc's
       rij2 += dr[k]*dr[k];
    }
    // compute inter-particle forces within cut-off distance:
    if (rij2 <= rcut2) {
       double rij = sqrt(rij2);             // center-to-center dist.
       double s2s = rij - (rad[i]+rad[j]);  // surface-to-surface dist.
       double fij = 2.5*K*pow((rcut-rij),1.5);    // Hertz contact force
       for (int k=0; k<3; k++) {
          f[i*3+k] += fij*dr[k]/rij;
          f[j*3+k] -= fij*dr[k]/rij;
       }
    }
}
