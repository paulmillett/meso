
# include "FIPISystem.hpp"
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

FIPISystem::FIPISystem(const Params& pin) : p(pin), c(p), k1(p), k2(p), k4(p),
                                            u(p), particles(p)
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
   // set some parameters:
   //	---------------------------------------

   visc = p.visc;    // fluid viscosity

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

FIPISystem::~FIPISystem()
{

}



// -------------------------------------------------------------------------
// Initialize the FIPI system:
// -------------------------------------------------------------------------

void FIPISystem::initializeFIPI()
{

   //	---------------------------------------
   // initialize the concentration field:
   //	---------------------------------------

   srand(time(NULL)*(p.rank+1));   // set the random seed
   double co = p.co;
   for (int i=0; i<nxyz; i++) {
      double r = (double)rand()/RAND_MAX;
      double val = co + 0.1*(r-0.5);
      c.setValue(i,val);
   }

   //	---------------------------------------
   // initialize the k1, k2, k4 fields:
   //	---------------------------------------

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
 				double kx = ii*(2.0*M_PI/p.Lx);
 				double ky = jj*(2.0*M_PI/p.Ly);
 				double kz = kk*(2.0*M_PI/p.Lz);
 				if (ii > p.NX/2) kx = (ii-p.NX)*(2.0*M_PI/p.Lx);
 				if (jj > p.NY/2) ky = (jj-p.NY)*(2.0*M_PI/p.Ly);
            if (kk > p.NZ/2) kz = (kk-p.NZ)*(2.0*M_PI/p.Lz);
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

   //	---------------------------------------
   // initialize the velocity field:
   //	---------------------------------------

   // something here...

   //	---------------------------------------
   // initialize the particle suspension:
   //	---------------------------------------

   particles.initParticles();

   //	---------------------------------------
   // Output the initial configuration:
   //	---------------------------------------

   int step = 0;
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   c.writeVTKFile("c",step,iskip,jskip,kskip);
   particles.writeVTKFile("part",step);

}



// -------------------------------------------------------------------------
// Take one step forward in the FIPI simulation:
// -------------------------------------------------------------------------

void FIPISystem::updateFIPI()
{

   double w = p.w;
   double kap = p.kap;
   particles.setTimeStep(current_step);


   // --------------- My Free-Energy Method below ----------------------

   // //	---------------------------------------
   // // update particles:
   // //	---------------------------------------
   //
   // Sfield eta = particles.mapToGrid();
   // if (current_step%10 == 0) {
   //    particles.zeroParticleForces();
   //    particles.calcCapillaryForce(eta,c);
   //    particles.pairwiseForces();
   //    particles.moveParticles();
   // }
   //
   // //	---------------------------------------
   // // update Cahn-Hilliard:
   // //	---------------------------------------
   //
   // Sfield dfdc = w*(c*c*c - c) + 2.0*eta*c;
   // c.fft(p_forward);
   // dfdc.fft(p_forward);
   // c = (c - p.dt*k2*dfdc)/(1.0 + kap*p.dt*k4);
   // c.ifft(p_backward);
   //
   // //	---------------------------------------
   // // calculate dfdc & grad_c:
   // //	---------------------------------------
   //
   // //Sfield dfdc = w*(1.0 - eta*(0.99)) * (c*c*c - c);
   // //Vfield grad_c = I*k1*c;
   // //grad_c.ifft(p_backward);
   // //c.ifft(p_backward);
   // //Vfield fip = particles.calcParticleInterfaceForce(c,grad_c);





   // --------------- My Pinning Force Method below ----------------------




   //	---------------------------------------
   // calculate dfdc & grad_c:
   //	---------------------------------------

   Sfield dfdc = w*(c*c*c - c);
   c.fft(p_forward);
   dfdc.fft(p_forward);
   Vfield grad_c = I*k1*c;
   grad_c.ifft(p_backward);
   c.ifft(p_backward);

   //	---------------------------------------
   // update particles:
   //	---------------------------------------

   particles.zeroParticleForces();
   Vfield fip = particles.calcParticleInterfaceForce(c,grad_c);
   particles.pairwiseForces();
   if (current_step%5 == 0) {
      particles.pairwiseForces();
      particles.moveParticles();
   }

   //	---------------------------------------
   // update Cahn-Hilliard:
   //	---------------------------------------

   Vfield fipc = fip*c;
   fipc.fft(p_forward);
   Sfield conv_c = fipc.div(I*k1);

   c.fft(p_forward);
   c = (c - p.dt*k2*dfdc - p.dt*conv_c)/(1.0 + kap*p.dt*k4);
   c.ifft(p_backward);




   // --------------- Botto's Method below ----------------------




   // //	---------------------------------------
   // // calculate dfdc, grad_c, and mu:
   // // [note: 'k2' does not have i^2, so
   // //  laplacians are multiplied by -1]
   // //	---------------------------------------
   //
   // Sfield dfdc = w*(c*c*c - c);
   // c.fft(p_forward);
   // dfdc.fft(p_forward);
   // Vfield grad_c = I*k1*c;
   // grad_c.ifft(p_backward);
   // Sfield conv_c = u.dot(grad_c);
   // conv_c.fft(p_forward);
   // //Sfield mu = dfdc + kap*k2*c;
   //
   // //	---------------------------------------
   // // update Cahn-Hilliard field:
   // //	---------------------------------------
   //
   // c = (c - p.dt*k2*dfdc - p.dt*conv_c)/(1.0 + kap*p.dt*k4);
   // c.ifft(p_backward);
   //
   // //	---------------------------------------
   // // update particles:
   // //	---------------------------------------
   //
   // particles.zeroParticleForces();
   // Vfield f_total = particles.calcParticleInterfaceForce(c,grad_c);
   // if (current_step < 4000 || current_step%10 == 0) {
   //    particles.pairwiseForces();
   //    particles.moveParticles();
   // }
   //
   // //	---------------------------------------
   // // calculate velocity field:
   // //	---------------------------------------
   //
   // // mu.ifft(p_backward);
   // // f_total += mu*grad_c;
   //
   // f_total.fft(p_forward);
   // u.resetVfield("f");
   // u = f_total - k1*(k1.dot(f_total))/k2;
   // u = u/(visc*k2);
   // u.ifft(p_backward);

}



// -------------------------------------------------------------------------
// Write output files:
// -------------------------------------------------------------------------

void FIPISystem::writeOutputFiles(int step)
{
   int iskip = p.iskip;
   int jskip = p.jskip;
   int kskip = p.kskip;
   c.writeVTKFile("c",step,iskip,jskip,kskip);
   particles.writeVTKFile("part",step);
   //u.writeVTKFile("u",step,iskip,jskip,kskip);
}
