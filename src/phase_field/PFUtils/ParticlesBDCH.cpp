
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
    // vector dimensions:
    for (int i=0; i<3*N; i++) {
        fcap.push_back(0.0);
        fcapSum.push_back(0.0);
    }
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

    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            double rr = (double)rand()/RAND_MAX;
            f[i*3+j] += sqrt(drag_coef*bm_str*2.0)*2.0*(rr-0.5);
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
    // MASTER broadcasts particle positions to everyone:
    MPI::COMM_WORLD.Bcast(&r[0],r.size(),MPI::DOUBLE,0);

    // loop over particles:
    Sfield eta(p);
    for (int pp=0; pp<N; pp++) {
        // width of region surrounding particles to probe:
        int wdth = int(ceil(2*rad[pp]+1));
        // particle's position:
        double x = r[pp*3+0];
        double y = r[pp*3+1];
        double z = r[pp*3+2];
        // decide if particle is in my domain:
        if (isMyParticle(x,pp) == true) {
            // grid point nearest particle (rounded down):
            int x0 = int(floor(x)/p.dx);
            int y0 = int(floor(y)/p.dy);
            int z0 = int(floor(z)/p.dz);
            // loop over region near particle:
            for (int i=0; i<wdth; i++) {
                int ii = x0 - (wdth/2 - 1) + i;
                if (ii < 0) ii += p.NX;
                if (ii > p.NX-1) ii -= p.NX;
                int ii_local = ii-p.xOff;
                if (inMyDomain(ii) == true) {
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
                            double val = 1.0 - r2/(rad[pp]*rad[pp]);
                            if (val < 0.0) val = 0.0;
                            double val0 = creal(eta.getValue(ii_local*p.nz*p.ny + jj*p.nz + kk));
                            eta.setValue(ii_local*p.nz*p.ny + jj*p.nz + kk, max(val,val0));
                        }
                    }
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

    // loop over particles:
    for (int pp=0; pp<N; pp++) {
        // width of region surrounding particles to probe:
        int wdth = int(ceil(2*rad[pp]+1));
        // particle's position:
        double x = r[pp*3+0];
        double y = r[pp*3+1];
        double z = r[pp*3+2];
        // decide if particle is in my domain:
        if (isMyParticle(x,pp) == true) {
            // grid point nearest particle (rounded down):
            int x0 = int(floor(x)/p.dx);
            int y0 = int(floor(y)/p.dy);
            int z0 = int(floor(z)/p.dz);
            // loop over region near particle:
            for (int i=0; i<wdth; i++) {
                int ii = x0 - (wdth/2 - 1) + i;
                if (ii < 0) ii += p.NX;
                if (ii > p.NX-1) ii -= p.NX;
                int ii_local = ii - p.xOff;
                if (inMyDomain(ii) == true) {
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
                            int ndx = ii_local*p.nz*p.ny + jj*p.nz + kk;
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
                            // should val = ccp?????
                            double val = 1.0 - rr*rr/(rad[pp]*rad[pp]);
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
    }
    // Collect all capillary forces to MASTER
    MPI::COMM_WORLD.Reduce(&fcap[0],&fcapSum[0],3*N,MPI::DOUBLE,MPI::SUM,0);
    if (p.rank == 0) {
        for (int i=0; i<3*N; i++) fcap[i] = fcapSum[i];
    }
    // make sure all ranks finish before moving on
    MPI::COMM_WORLD.Barrier();
}



// -------------------------------------------------------------------------
// Function to determine if particle's x-position is inside
// my domain (or within a certain halo distance):
// -------------------------------------------------------------------------

bool ParticlesBDCH::isMyParticle(double x, int i)
{
    double lower = double(p.xOff)*p.dx - rad[i];
    double upper = double(p.xOff)*p.dx + double(p.nx)*p.dx + rad[i];
    if (x >= lower && x <= upper)
        return true;
    // implement PBC for first and last domain
    else if ( p.rank == 0 && x > (lower+double(p.NX)*p.dx) )
        return true;
    else if ( p.rank == p.np && x < rad[i] )
        return true;
    else
        return false;
}



// -------------------------------------------------------------------------
// Function to determine if grid x-position is inside
// my domain:
// -------------------------------------------------------------------------

bool ParticlesBDCH::inMyDomain(int i)
{
    if (i >= p.xOff && i <= p.xOff+p.nx-1)
        return true;
    else
        return false;
}
