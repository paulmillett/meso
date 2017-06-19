
# include "ParticlesBDCH.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

ParticlesBDCH::ParticlesBDCH(const CommonParams& pin,
                             const GetPot& input_params) :
                             PDParticles(pin,input_params), p(pin)
{
    // get the capillary force parameters:
    cap_str = input_params("PFApp/cap_str",1.0);
    cout << cap_str << endl;
    // vector dimensions:
    for (int i=0; i<N; i++)
        for (int k=0; k<3; k++)
            fcap.push_back(0.0);
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

ParticlesBDCH::~ParticlesBDCH()
{

}



// -------------------------------------------------------------------------
// Auxiliary forces:
// -------------------------------------------------------------------------

void ParticlesBDCH::auxiliaryForces()
{

    cout << "Hi" << endl;

    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            double rr = (double)rand()/RAND_MAX;
            f[i*3+j] += bm_str*2.0*(rr-0.5);
            f[i*3+j] -= drag_coef*v[i*3+j];
            f[i*3+j] += fcap[i*3+j];
        }
    }
}



// -------------------------------------------------------------------------
// Map particle positions onto an Sfield grid:
// -------------------------------------------------------------------------

Sfield ParticlesBDCH::mapToGrid()
{
    // width of region surrounding particles to probe:
    int wdth = 7;
    // loop over particles:
    Sfield eta(p);
    for (int pp=0; pp<N; pp++) {
        // particle's position:
        double x = r[pp*3+0];
        double y = r[pp*3+1];
        double z = r[pp*3+2];
        // grid point nearest particle (rounded down):
        int x0 = int(floor(x)/p.dx);
        int y0 = int(floor(y)/p.dy);
        int z0 = int(floor(z)/p.dz);
        // loop over region near particle:
        for (int i=0; i<wdth; i++) {
            int ii = x0 - (wdth/2 - 1) + i;
            if (ii < 0) ii += p.NX;
            if (ii > p.NX-1) ii -= p.NX;
            for (int j=0; j<wdth; j++) {
                int jj = y0 - (wdth/2 - 1) + j;
                if (jj < 0) jj += p.NY;
                if (jj > p.NY-1) jj -= p.NY;
                for (int k=0; k<wdth; k++) {
                    int kk = z0 - (wdth/2 - 1) + k;
                    if (kk < 0) kk += p.NZ;
                    if (kk > p.NZ-1) kk -= p.NZ;
                    if (p.NZ == 1) kk = 0;
                    // calculate distance to point:
                    double rx = x - double(ii);
                    double ry = y - double(jj);
                    double rz = z - double(kk);
                    rx -= round(rx/(p.NX*p.dx))*p.NX*p.dx;
                    ry -= round(ry/(p.NY*p.dy))*p.NY*p.dy;
                    rz -= round(rz/(p.NZ*p.dz))*p.NZ*p.dz;
                    double r2 = rx*rx + ry*ry + rz*rz;
                    // assign spread function to grid:
                    double val = 1.0 - r2/6.25;
                    if (val < 0.0) val = 0.0;
                    double val0 = creal(eta.getValue(ii*p.nz*p.ny + jj*p.nz + kk));
                    eta.setValue(ii*p.nz*p.ny + jj*p.nz + kk, max(val,val0));
                }
            }
        }
    }
    return eta;
}



// -------------------------------------------------------------------------
// Calculate capillary forces on particles due to fluid-fluid interface:
// -------------------------------------------------------------------------

void ParticlesBDCH::calcCapillaryForce(const Sfield& cp,
                                       const Sfield& c1,
                                       const Sfield& c2)
{
    // zero fcap array:
    for (int i=0; i<3*N; i++) fcap[i] = 0.0;
    // width of region surrounding particles to probe:
    int wdth = 7;
    // loop over particles:
    for (int pp=0; pp<N; pp++) {
        // particle's position:
        double x = r[pp*3+0];
        double y = r[pp*3+1];
        double z = r[pp*3+2];
        // grid point nearest particle (rounded down):
        int x0 = int(floor(x)/p.dx);
        int y0 = int(floor(y)/p.dy);
        int z0 = int(floor(z)/p.dz);
        // loop over region near particle:
        for (int i=0; i<wdth; i++) {
            int ii = x0 - (wdth/2 - 1) + i;
            if (ii < 0) ii += p.NX;
            if (ii > p.NX-1) ii -= p.NX;
            for (int j=0; j<wdth; j++) {
                int jj = y0 - (wdth/2 - 1) + j;
                if (jj < 0) jj += p.NY;
                if (jj > p.NY-1) jj -= p.NY;
                for (int k=0; k<wdth; k++) {
                    int kk = z0 - (wdth/2 - 1) + k;
                    if (kk < 0) kk += p.NZ;
                    if (kk > p.NZ-1) kk -= p.NZ;
                    if (p.NZ == 1) kk = 0;
                    // get local field values:
                    int ndx = ii*p.nz*p.ny + jj*p.nz + kk;
                    double ccp = creal(cp.getValue(ndx));
                    double cc1 = creal(c1.getValue(ndx));
                    double cc2 = creal(c2.getValue(ndx));
                    // calculate distance to point:
                    double rx = double(ii) - x;
                    double ry = double(jj) - y;
                    double rz = double(kk) - z;
                    rx -= round(rx/(p.NX*p.dx))*p.NX*p.dx;
                    ry -= round(ry/(p.NY*p.dy))*p.NY*p.dy;
                    rz -= round(rz/(p.NZ*p.dz))*p.NZ*p.dz;
                    double rr = sqrt(rx*rx + ry*ry + rz*rz);
                    // calculate interface force on particle:
                    double val = 1.0 - rr*rr/6.25;
                    if (val < 0.0) val = 0.0;
                    double fint = val*cc1*cc2;
                    fcap[pp*3+0] += cap_str*fint*(rx);///rr);
                    fcap[pp*3+1] += cap_str*fint*(ry);///rr);
                    fcap[pp*3+2] += cap_str*fint*(rz);///rr);
                }
            }
        }
    }
}
