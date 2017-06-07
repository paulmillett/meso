
# include "Hertz.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

Hertz::Hertz(const CommonParams& pin,
             const GetPot& input_params) : PDBaseClass(pin,input_params),
                                           p(pin)
{

   //	---------------------------------------
   // set needed parameters:
   //	---------------------------------------

   K = input_params("PDApp/K",0.5);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

Hertz::~Hertz()
{

}



// -------------------------------------------------------------------------
// Function that calculates pairwise force based on separation distance:
// -------------------------------------------------------------------------

void Hertz::fijFunc(int i, int j)
{
    // calculate separation vector rij of particle i and j
    for (size_t k=0;k<3;k++) rij[k] = r[i*3+k]-r[j*3+k];
    // calculate its magnitude
    rijMag = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
    // calculate its unit vector
    for (size_t k=0;k<3;k++) rijUnit[k] = rij[k]/rijMag;
    // calculate the surface-to-surface distance
    s2s = rijMag - rad[i] - rad[j];
    
    // calculate the magnitude of the force
    if (s2s < 0.0) fmag = K*pow(-s2s,1.5);
    else fmag = 0.0;

    // assign forces to particles i and j
    for (size_t k=0;k<3;k++) f[i*3+k] = fmag*rijUnit[k];
    for (size_t k=0;k<3;k++) f[j*3+k] = -fmag*rijUnit[k];
}
