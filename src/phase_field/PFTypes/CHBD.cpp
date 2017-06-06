
# include "CHBD.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHBD::CHBD(const CommonParams& pin,
           const GetPot& input_params) : p(pin), c1(p), c2(p), particles(p),
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

CHBD::~CHBD()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void CHBD::initPhaseField()
{

   //	---------------------------------------
   // initialize the concentration field:
   //	---------------------------------------

   srand(time(NULL)*(p.rank+1));   // set the random seed
   for (int i=0; i<nxyz; i++) {
      double r = (double)rand()/RAND_MAX;
      double val = co + 0.1*(r-0.5);
      c1.setValue(i,val);
      c2.setValue(i,1-val);
   }

   //	---------------------------------------
   // initialize the particle suspension:
   //	---------------------------------------

   particles.initParticles();

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

void CHBD::updatePhaseField()
{
   //	---------------------------------------
   // Update particles:
   //	---------------------------------------

   Sfield cp = particles.mapToGrid();
   if (current_step%20 == 0) {
      particles.zeroParticleForces();
      particles.calcCapillaryForce(cp,c1,c2);
      particles.pairwiseForces();
      particles.moveParticles();
   }

   //	---------------------------------------
   // Update Cahn-Hilliard:
   //	---------------------------------------

   Sfield dfdc1 = 12.0*w*(c1*c1*c1 - c1*c1 + c1*c2*c2 + c1*cp);
   Sfield dfdc2 = 12.0*w*(c2*c2*c2 - c2*c2 + c2*c1*c1 + c2*cp);
   c1.fft(p_forward);
   c2.fft(p_forward);
   dfdc1.fft(p_forward);
   dfdc2.fft(p_forward);
   c1 = (c1 - p.dt*k2*dfdc1)/(1.0 + kap*p.dt*k4);
   c2 = (c2 - p.dt*k2*dfdc2)/(1.0 + kap*p.dt*k4);
   c1.ifft(p_backward);
   c2.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void CHBD::outputPhaseField()
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   Sfield phi = c1 - c2;
   phi.writeVTKFile("c",current_step,iskip,jskip,kskip);
   particles.writeVTKFile("part",current_step);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void CHBD::calculateKfields()
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
