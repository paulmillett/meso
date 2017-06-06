
# include "TIPS2.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPS2::TIPS2(const CommonParams& pin,
             const GetPot& input_params) : p(pin), c(p), phi(p),
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

   co = input_params("PFApp/co",0.3);
   w = input_params("PFApp/w",1.0);
   M = input_params("PFApp/M",1.0);
   kapc = input_params("PFApp/kapc",1.0);
   kapp = input_params("PFApp/kapp",1.0);
   phio = input_params("PFApp/phio",0.5);
   prate = input_params("PFApp/prate",0.00001);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPS2::~TIPS2()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPS2::initPhaseField()
{

   //	---------------------------------------
   // initialize the concentration field:
   //	---------------------------------------

   srand(time(NULL)*(p.rank+1));   // set the random seed
   for (int i=0; i<nxyz; i++) {
      double r1 = (double)rand()/RAND_MAX;
      double r2 = (double)rand()/RAND_MAX;
      double val1 = co + 0.4*(r1-0.5);
      double val2 = phio + 0.1*(r2-0.5);
      c.setValue(i,val1);
      phi.setValue(i,val2);
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

void TIPS2::updatePhaseField()
{
   // calculate positions on phase diagram:
   double t = double(current_step);
   double cA = 0.1; //0.1 - (t/10000.0)*0.1;
   double cB = 0.55; //0.4 + (t/10000.0)*0.4;
   // update 'phi':
   Sfield B = -5.0*(c-0.1)*(c-0.1)*(c-0.1);
   Sfield dphi = w*(phi*phi*phi - phi + B*(1 - 2*phi*phi + phi*phi*phi*phi));
   phi.fft(p_forward);
   dphi.fft(p_forward);
   phi = (phi - p.dt*k2*dphi)/(1.0 + kapp*p.dt*k4);
   phi.ifft(p_backward);
   phi -= p.dt*prate;
   // update 'c':
   Sfield A = 0.5*(phi + 1);
   Sfield D = 1.0*(1-A) + 0.1*A;
   Sfield dfdc = 2*D*c;
   c.fft(p_forward);
   dfdc.fft(p_forward);
   c = (c - M*p.dt*k2*dfdc)/(1.0 + M*kapc*p.dt*k4);
   c.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPS2::outputPhaseField()
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   c.writeVTKFile("c",current_step,iskip,jskip,kskip);
   phi.writeVTKFile("phi",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void TIPS2::calculateKfields()
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
