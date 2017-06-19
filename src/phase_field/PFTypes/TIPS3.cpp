
# include "TIPS3.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPS3::TIPS3(const CommonParams& pin,
             const GetPot& input_params) : p(pin), c(p), k1(p), k2(p), k4(p)
{

    //	---------------------------------------
    // create the FFTw plans:
    //	---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    dummy = fftw_alloc_complex(nxyz);
    p_forward  = fftw_mpi_plan_dft_3d(p.NX,p.NY,p.NZ,dummy,dummy,
                 MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE);
    p_backward = fftw_mpi_plan_dft_3d(p.NX,p.NY,p.NZ,dummy,dummy,
                 MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE);
    fftw_free(dummy);

    //	---------------------------------------
    // set needed parameters:
    //	---------------------------------------

    co = input_params("PFApp/co",0.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    chi = input_params("PFApp/chi",1.0);
    kT = input_params("PFApp/kT",1.0);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    B = input_params("PFApp/A",0.4);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPS3::~TIPS3()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPS3::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=0; i<nxyz; i++) {
        double r = (double)rand()/RAND_MAX;
        double val = co + 0.1*(r-0.5);
        c.setValue(i,val);
    }

    //	---------------------------------------
    // initialize the fourier wave-vectors:
    //	---------------------------------------

    calculateKfields();

    //	---------------------------------------
    // Output the initial configuration:
    //	---------------------------------------

    current_step = 0;
    outputPhaseField();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void TIPS3::updatePhaseField()
{
    // calculate energy parameters:
    double scl = 7.0 - (7.0-4.0)/double(p.nstep)*double(current_step);
    // calculate chemical potential from Flory-Huggins theory:
    Sfield dfdc(p);
    for (int i=0; i<p.nx; i++) {
        for (int j=0; j<p.ny; j++) {
            for (int k=0; k<p.nz; k++) {
                int ndx = i*p.nz*p.ny + j*p.nz + k;
                double cc = creal(c.getValue(ndx));
                //  double df = (log(cc) + 1.0)/N - log(1.0-cc) - 1.0 +
                //              chi*(1.0-2.0*cc);
                //  df *= kT;
                //  if (cc <= 0.0) df = -1.5*A*sqrt(-cc);
                double df = 4*scl*cc*cc*cc - 12*cc*cc;
                if (cc <= 0.0) df = -1.5*A*sqrt(-cc);
                dfdc.setValue(ndx,df);
            }
        }
    }
    // update c:
    c.fft(p_forward);
    dfdc.fft(p_forward);
    c = (c - M*p.dt*k2*dfdc)/(1.0 + M*kap*p.dt*k4);
    c.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPS3::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void TIPS3::calculateKfields()
{
    Sfield kxf(p);
    Sfield kyf(p);
    Sfield kzf(p);
    for (int i=0; i<p.nx; i++) {
        for (int j=0; j<p.ny; j++) {
            for (int k=0; k<p.nz; k++) {
                int ii = i + p.xOff;
                int jj = j;
                int kk = k;
                int ndx = i*p.NY*p.NZ + j*p.NZ + k;
                double kx = ii*(2.0*M_PI/p.LX);
                double ky = jj*(2.0*M_PI/p.LY);
                double kz = kk*(2.0*M_PI/p.LZ);
                if (ii > p.NX/2) kx = (ii-p.NX)*(2.0*M_PI/p.LX);
                if (jj > p.NY/2) ky = (jj-p.NY)*(2.0*M_PI/p.LY);
                if (kk > p.NZ/2) kz = (kk-p.NZ)*(2.0*M_PI/p.LZ);
                double kx2 = kx*kx;
                double ky2 = ky*ky;
                double kz2 = kz*kz;
                double kijk2 = kx2 + ky2 + kz2;
                double kijk4 = kijk2*kijk2;
                kxf.setValue(ndx,kx);
                kyf.setValue(ndx,ky);
                kzf.setValue(ndx,kz);
                k2.setValue(ndx,kijk2);
                k4.setValue(ndx,kijk4);
            }
        }
    }
    k1.setXValues(kxf);
    k1.setYValues(kyf);
    k1.setZValues(kzf);
}
