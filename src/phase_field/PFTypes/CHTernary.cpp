
# include "CHTernary.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHTernary::CHTernary(const CommonParams& pin,
                     const GetPot& input_params) : p(pin), c1(p), c2(p), c3(p),
                                                   k1(p), k2(p), k4(p)
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
   w = input_params("PFApp/w",1.0);
   M = input_params("PFApp/M",1.0);
   kap = input_params("PFApp/kap",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHTernary::~CHTernary()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void CHTernary::initPhaseField()
{

   //	---------------------------------------
   // initialize the concentration field:
   //	---------------------------------------

   srand(time(NULL)*(p.rank+1));   // set the random seed
   for (int i=0; i<nxyz; i++) {
      double r1 = (double)rand()/RAND_MAX;
      double r2 = (double)rand()/RAND_MAX;
      double r3 = (double)rand()/RAND_MAX;
      double val1 = 0.47 + 0.1*(r1-0.5);
      double val2 = 0.47 + 0.1*(r2-0.5);
      double val3 = 0.10 + 0.1*(r3-0.5);
      c1.setValue(i,val1);
      c2.setValue(i,val2);
      c3.setValue(i,val3);
   }

   //	---------------------------------------
   // initialize the fourier wave-vectors:
   //	---------------------------------------

   calculateKfields();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void CHTernary::updatePhaseField()
{
   //Sfield dfdc1 = w*(4*c1*c1*c1 - 6*c1*c1 + 2*c1 + 2*c1*(c2*c2 + c3*c3));
   //Sfield dfdc2 = w*(4*c2*c2*c2 - 6*c2*c2 + 2*c2 + 2*c2*(c1*c1 + c3*c3));
   //Sfield dfdc3 = w*(4*c3*c3*c3 - 6*c3*c3 + 2*c3 + 2*c3*(c1*c1 + c2*c2));
   Sfield dfdc1 = w*(c1*c1*c1 - c1 + 2*c1*(c2*c2 + 0.5*c3*c3));
   Sfield dfdc2 = w*(c2*c2*c2 - c2 + 2*c2*(c1*c1 + 0.5*c3*c3));
   // Sfield dfdc3 = w*(c3*c3*c3 - c3 + 2*c3*(c1*c1 + c2*c2)*0.5);
   // for (int i=0; i<nxyz; i++) {
   //    double c3val = creal(c3.getValue(i));
   //    if (c3val < 0.0) dfdc3.addValue(i,w*2*c3val);
   // }
   Sfield dfdc3 = (4*c3*c3*c3 - 6*c3*c3 + 2*c3)*(1-3*c1*c2);
   c1.fft(p_forward);
   c2.fft(p_forward);
   c3.fft(p_forward);
   dfdc1.fft(p_forward);
   dfdc2.fft(p_forward);
   dfdc3.fft(p_forward);
   c1 = (c1 - p.dt*k2*dfdc1)/(1.0 + kap*p.dt*k4);
   c2 = (c2 - p.dt*k2*dfdc2)/(1.0 + kap*p.dt*k4);
   c3 = (c3 - p.dt*k2*dfdc3)/(1.0 + kap*p.dt*k4);
   c1.ifft(p_backward);
   c2.ifft(p_backward);
   c3.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void CHTernary::outputPhaseField()
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   Sfield c12 = c1 - c2;
   c12.writeVTKFile("c12",current_step,iskip,jskip,kskip);
   c3.writeVTKFile("c3",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void CHTernary::calculateKfields()
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
