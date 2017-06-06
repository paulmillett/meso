
# include "CHSystem.hpp"
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHSystem::CHSystem(const Params& pin) : p(pin)
{

	//	---------------------------------------
   //	Index offset for stencil neighbors:
   //	---------------------------------------

	deli = (p.nz+2)*(p.ny+2);
	delj = (p.nz+2);
	delk = 1;

   //	---------------------------------------
   //	Create conc. & chem. potential vectors:
   //	---------------------------------------

	for (int i=0; i<p.nc; i++) {
		c.push_back(FDVector(p,0.0));
		mu.push_back(FDVector(p,0.0));
	}

   //	---------------------------------------
   //	Create objects for various purposes:
   //	---------------------------------------

   muObj = ChemPot::ChemPotFactory(p,c);
	srcObj = SourceTerms::SourceTermsFactory(p,c);
	mobObj = Mobilities::MobilitiesFactory(p,c);

	//	---------------------------------------
	//	Get some extra data:
	//	---------------------------------------

 	mu_noise = p.input_params("CHApp/free_energy/noise",0.0);
	src_type = p.input_params("CHApp/source/type","none");

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHSystem::~CHSystem()
{
   delete muObj;
	delete srcObj;
	delete mobObj;
}



// -------------------------------------------------------------------------
// Initialize the Cahn-Hilliard system:
// -------------------------------------------------------------------------

void CHSystem::initializeCahnHilliard()
{
	double c0_init_mean = p.input_params("CHApp/initial_condition/c0_init_mean",0.5);
	double c0_init_noise = p.input_params("CHApp/initial_condition/c0_init_noise",0.1);
	c[0].initMeanRandom(c0_init_mean,c0_init_noise);
}



// -------------------------------------------------------------------------
// Take one step forward in the Cahn-Hilliard simulation:
// -------------------------------------------------------------------------

void CHSystem::updateCahnHilliard()
{

   // ----------------------------------------
   //	Periodic boundary conditions for 'c':
   // ----------------------------------------

	for (int i=0; i<p.nc; i++) {
		c[i].updateBoundaryConditions();
	}
	MPI::COMM_WORLD.Barrier();

   // ----------------------------------------
   //	Calculate chemical potential:
   // ----------------------------------------

	for (int n=0; n<p.nc; n++) {
		for (int i=1; i<p.nx+1; i++) {
			for (int j=1; j<p.ny+1; j++) {
				for (int k=1; k<p.nz+1; k++) {
					int ndx = k*delk + j*delj + i*deli;
					double chempot = muObj->muFunc(ndx,n,i,j,k);
         		mu[n].setValue(ndx,chempot);
         	}
         }
		}
	}

	// ----------------------------------------
   //	Add random thermal fluctuations to 'mu':
   // ----------------------------------------

	if (mu_noise > 0.0) {
		for (int n=0; n<p.nc; n++) {
			for (int i=1; i<p.nx+1; i++) {
				for (int j=1; j<p.ny+1; j++) {
					for (int k=1; k<p.nz+1; k++) {
						int ndx = k*delk + j*delj + i*deli;
						double r = (double)rand()/RAND_MAX;
						mu[n].addValue(ndx,mu_noise*(0.5-r));
					}
				}
			}
		}
	}

   // ----------------------------------------
   //	Periodic boundary conditions for 'mu':
   // ----------------------------------------

	for (int i=0; i<p.nc; i++) {
		mu[i].updateBoundaryConditions();
	}
	MPI::COMM_WORLD.Barrier();

   // ----------------------------------------
   //	Update C-H equation:
   // ----------------------------------------

	for (int n=0; n<p.nc; n++) {
		for (int i=1; i<p.nx+1; i++) {
			for (int j=1; j<p.ny+1; j++) {
         	for (int k=1; k<p.nz+1; k++) {
					int ndx = k*delk + j*delj + i*deli;
					double chRHS = mu[n].calculateLaplacian(ndx);
            	c[n].updateExplicitEuler(ndx,p.dt,chRHS);
            }
         }
		}
	}

	// ----------------------------------------
   //	Add extra 'source' terms to solution:
   // ----------------------------------------

	if (src_type != "none") {
		for (int n=0; n<p.nc; n++) {
			for (int i=1; i<p.nx+1; i++) {
				for (int j=1; j<p.ny+1; j++) {
					for (int k=1; k<p.nz+1; k++) {
						int ndx = k*delk + j*delj + i*deli;
						double source = srcObj->srcFunc(ndx,n,i,j,k);
						c[n].updateExplicitEuler(ndx,p.dt,source);
					}
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Print VTK output file:
// -------------------------------------------------------------------------

void CHSystem::writeOutputFiles(int step)
{
	int iskip = p.input_params("Output/iskip",1);
   int jskip = p.input_params("Output/jskip",1);
   int kskip = p.input_params("Output/kskip",1);
	c[0].writeVTKFile("c",step,iskip,jskip,kskip);
}
