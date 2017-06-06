
# include "CHFFTSystem.hpp"
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHFFTSystem::CHFFTSystem(const Params& pin) : p(pin)
{

	//	---------------------------------------
	//	Index offset for stencil neighbors:
	//	---------------------------------------

	deli = (2*(p.nz/2+1))*p.ny;
	delj = (2*(p.nz/2+1));
	delk = 1;

	//	---------------------------------------
	//	Get local array size:
	//	---------------------------------------

	sizeR = p.nx*p.ny*(2*(p.nz/2+1));

   //	---------------------------------------
   //	Create array of chemical species:
   //	---------------------------------------

	c = new FFTArrays*[p.nc];
   for (int i=0; i<p.nc; i++) {
      c[i] = new FFTArrays(p);
   }

	//	---------------------------------------
	//	Create auxillary vector 'eta':
	//	---------------------------------------

	for (int i=0; i<sizeR; i++) {
		eta.push_back(0.0);
	}

	//	---------------------------------------
   //	Create objects for various purposes:
   //	---------------------------------------

	icObj = FFTInitCond::FFTInitCondFactory(p,c,eta);
	muObj = FFTChemPot::FFTChemPotFactory(p,c,eta);
	srcObj = FFTSourceTerms::FFTSourceTermsFactory(p,c);
	mobObj = FFTMobilities::FFTMobilitiesFactory(p,c);

	//	---------------------------------------
	//	Get some extra data:
	//	---------------------------------------

	c_noise  = p.input_params("CHFFTApp/source/noise",0.0);
	mu_noise = p.input_params("CHFFTApp/free_energy/noise",0.0);
	src_type = p.input_params("CHFFTApp/source/type","none");
	mob_type = p.input_params("CHFFTApp/mobility/type","isotropic");

	//	---------------------------------------
	//	Create boolian array that indicates
	// whether mobility is isotropic:
	//	---------------------------------------

	mob_flag = new bool[p.nc];
	for (int i=0; i<p.nc; i++) {
		if (mob_type == "isotropic") mob_flag[i] = true;
		if (mob_type != "isotropic") mob_flag[i] = false;
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHFFTSystem::~CHFFTSystem()
{
   for (int i=0; i<p.nc; i++) {      // deallocate 'c' array
      delete c[i];
   }
   delete[] c;
	delete[] mob_flag;
}



// -------------------------------------------------------------------------
// Initialize the Cahn-Hilliard system:
// -------------------------------------------------------------------------

void CHFFTSystem::initializeCahnHilliard()
{
	icObj->icFunc();
}



// -------------------------------------------------------------------------
// Take one step forward in the Cahn-Hilliard simulation:
// -------------------------------------------------------------------------

void CHFFTSystem::updateCahnHilliard()
{

   // ----------------------------------------
   //	Calculate chemical potential:
   // ----------------------------------------

   for (int n=0; n<p.nc; n++) {
		for (int i=0; i<p.nx; i++) {
			for (int j=0; j<p.ny; j++) {
				for (int k=0; k<p.nz; k++) {
					int ndx = i*deli + j*delj + k*delk;
					double cpot = muObj->muFunc(n,ndx,i,j,k);
					c[n]->setMu(ndx,cpot);
				}
			}
		}
   }

	// ----------------------------------------
   //	Add random thermal fluctuations to 'mu':
   // ----------------------------------------

	if (mu_noise > 0.0) {
      for (int n=0; n<p.nc; n++) {
         for (int i=0; i<sizeR; i++) {
				double r = (double)rand()/RAND_MAX;
				c[n]->addtoMu(i,mu_noise*(0.5-r));
         }
      }
   }

	// ----------------------------------------
   //	Calculate mobility function (only if
	// mob_flag is false):
   // ----------------------------------------

	for (int n=0; n<p.nc; n++) {
		if (mob_flag[n] == false) {
			for (int i=0; i<p.nx; i++) {
				for (int j=0; j<p.ny; j++) {
					for (int k=0; k<p.nz; k++) {
						int ndx = i*deli + j*delj + k*delk;
						double mob = mobObj->mobFunc(n,ndx,i,j,k);
						c[n]->setMob(ndx,mob);
					}
				}
			}
		}
   }

	// ----------------------------------------
   //	Forward transform 'c' and 'mu':
   // ----------------------------------------

   for (int n=0; n<p.nc; n++) {
      c[n]->forwardFFT();
   }

	// ----------------------------------------
   //	Update 'c' in Fourier space:
   // ----------------------------------------

   for (int n=0; n<p.nc; n++) {
      c[n]->updateConcFFT(mob_flag[n]);
   }

   // ----------------------------------------
   //	Backward transform 'c' into real space:
   // ----------------------------------------

   for (int n=0; n<p.nc; n++) {
      c[n]->backwardFFT();
   }

	// ----------------------------------------
   //	Add extra 'source' terms to solution:
   // ----------------------------------------

	if (src_type != "none") {
		for (int n=0; n<p.nc; n++) {
         for (int i=0; i<sizeR; i++) {
				double src = srcObj->srcFunc(n,i);
				c[n]->addtoConc(i,p.dt*src);
         }
      }
	}

	// ----------------------------------------
   //	Add random fluctuations to 'c':
   // ----------------------------------------

	if (c_noise > 0.0) {
      for (int n=0; n<p.nc; n++) {
         for (int i=0; i<sizeR; i++) {
				double r = (double)rand()/RAND_MAX;
				c[n]->addtoConc(i,p.dt*c_noise*(0.5-r));
         }
      }
   }

}



// -------------------------------------------------------------------------
// Write output files:
// -------------------------------------------------------------------------

void CHFFTSystem::writeOutputFiles(int step)
{
	int iskip = p.input_params("Output/iskip",1);
   int jskip = p.input_params("Output/jskip",1);
   int kskip = p.input_params("Output/kskip",1);

	switch (p.nc) {
		case 1:
			c[0]->writeVTKFile("c",step,iskip,jskip,kskip);
			break;
		case 2:
			c[0]->writeVTKFile("c1",step,iskip,jskip,kskip);
			c[1]->writeVTKFile("c2",step,iskip,jskip,kskip);
			break;
	}
}
