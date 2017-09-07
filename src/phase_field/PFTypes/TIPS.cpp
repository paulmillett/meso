
# include "TIPS.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPS::TIPS(const CommonParams& pin,
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
   w = input_params("PFApp/w",1.0);
   M = input_params("PFApp/M",1.0);
   kap = input_params("PFApp/kap",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPS::~TIPS()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPS::initPhaseField()
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

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void TIPS::updatePhaseField()
{
   // calculate positions on phase diagram:
   double t = double(current_step);
   double cA = 0.1; //0.1 - (t/10000.0)*0.1;
   double cB = 0.55; //0.4 + (t/10000.0)*0.4;
   // calculate chemical potential:
   Sfield ccA = c - cA;
   Sfield ccB = c - cB;
   Sfield dfdc = 2.0*w*(ccA*ccB*ccB + ccA*ccA*ccB);
   // calculate mobility field:
   //Sfield cm = 1.0 - c;
   Sfield mob(p);// = cm*cm; //*cm*cm*cm*cm;
   for (int i=0; i<p.nx; i++) {
      for (int j=0; j<p.ny; j++) {
         for (int k=0; k<p.nz; k++) {
            //int ii = i + p.xOff;
            //double Mc = 0.5*(tanh(double(ii - 100)) + 1.0) +
            //            0.5*(1 - tanh(double(ii-30)));
            int ndx = i*p.nz*p.ny + j*p.nz + k;
            double cc = creal(c.getValue(ndx));
            double Mc = 1.0;
            if (cc > 0.1) Mc = 0.018/(pow(cc,1.75));
            mob.setValue(ndx,Mc);
         }
      }
   }

   // update concentration field:
   c.fft(p_forward);
   dfdc.fft(p_forward);
   Sfield mu = dfdc + kap*k2*c;
   Vfield gradmu = I*k1*mu;
   gradmu.ifft(p_backward);
   Vfield Mgradmu = mob*gradmu;
   Mgradmu.fft(p_forward);
   Sfield rhs = k1.dot(Mgradmu);
   rhs *= I;
   double A = 0.8;
   c = c + rhs*p.dt/(1.0 + A*p.dt*kap*k4);
   //c = (c - p.dt*k2*dfdc)/(1.0 + kap*p.dt*k4);
   c.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPS::outputPhaseField()
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void TIPS::calculateKfields()
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
