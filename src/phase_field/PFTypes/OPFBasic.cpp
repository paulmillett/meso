
# include "OPFBasic.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

OPFBasic::OPFBasic(const CommonParams& pin,
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
	S = input_params("PFApp/S",0.05);
	chiN = input_params("PFApp/chiN",15);
	chiN2 = chiN*chiN;
	chiN3 = chiN2*chiN;
	chiN4 = chiN2*chiN2;
   
	c_2 = S*(-0.00368*chiN2 - 1.964*chiN +15.99); 
	c_4 = S*(0.01882*chiN2 + 3.176*chiN -26.9); 
	c_5 = S*(-0.000003452*chiN4 + 0.0004019*chiN3 - 0.01791*chiN2 + 0.3735*chiN -1.684); 
	c_7 = S*(0.000004339*chiN4 - 0.000446*chiN3 + 0.01729*chiN2 - 0.32*chiN + 11.59); 
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

OPFBasic::~OPFBasic()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void OPFBasic::initPhaseField()
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
   //outputPhaseField();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void OPFBasic::updatePhaseField()
{
   Sfield dfdc = 2*c_2*(c-0.5) + 4*c_4*(c-0.5)*(c-0.5)*(c-0.5);
   c.fft(p_forward);
   dfdc.fft(p_forward);
   c = (c - p.dt*k2*dfdc)/(1.0 + c_5*p.dt*k4);
   c.ifft(p_backward);
   c -= p.dt*c_7*(c - co);
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void OPFBasic::outputPhaseField()
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate the wave-vector fields (k^1, k^2, k^4 for spectral solution):
// -------------------------------------------------------------------------

void OPFBasic::calculateKfields()
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
