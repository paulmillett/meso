
# include "scsp2D.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

scsp2D::scsp2D(const CommonParams& pin, const GetPot& input_params) : p(pin)
{

	// ---------------------------------------
	// set needed parameters:
	// ---------------------------------------

	tau = input_params("LBApp/tau",1.0);
	mu = input_params("LBApp/mu",1.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

scsp2D::~scsp2D()
{

}



// -------------------------------------------------------------------------
// Initialize lattice-boltzmann method:
// -------------------------------------------------------------------------

void scsp2D::initLatticeBoltzmann()
{

	// ---------------------------------------
	// initialize the concentration field:
	// ---------------------------------------

	srand(time(NULL)*(p.rank+1));   // set the random seed
	for (int i=0; i<nxyz; i++) {
		double r = (double)rand()/RAND_MAX;
		double val = co + 0.1*(r-0.5);
		c.setValue(i,val);
	}

	// ---------------------------------------
	// Output the initial configuration:
	// ---------------------------------------

	current_step = 0;
	outputLatticeBoltzmann();

}



// -------------------------------------------------------------------------
// Step forward in time the lattice-boltzmann method:
// -------------------------------------------------------------------------

void CHBasic::updateLatticeBoltzmann()
{

	// ---------------------------------------
	// update fprev = f:
	// ---------------------------------------

	for (int i=0; i<nxyzn; i++) fprev[i] = f[i];

	// ---------------------------------------
	// update macros:
	// ---------------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double sum  = 0.0;
				double sumx = 0.0;
				double sumy = 0.0;
				for (int n=0; n<nn; n++) {
					int ndxn = ndx*nn + n;
					sum += f[ndxn];
					sumx += f[ndxn]*ex[n];
					sumy += f[ndxn]*ey[n];
				}
				rho[ndx] = sum;
				u[ndx] = sumx/rho[ndx];
				v[ndx] = sumy/rho[ndx];
			}
		}
	}

	// ---------------------------------------
	// fluid forces:
	// ---------------------------------------


	// ---------------------------------------
	// equilibrium distributions:
	// ---------------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double u_F = u[ndx] + tau*fx[ndx]/rho[ndx];
				double v_F = v[ndx] + tau*fy[ndx]/rho[ndx];
				double uv2 = u_F*u_F + v_F*v_F;
				for (int n=0; n<nn; n++) {
					int ndxn = ndx*nn + n;
					double evel = ex[n]*u_F + ey[n]*v_F;
					feq[ndxn] = rho[ndx]*wa[n]*(1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uv2);
				}
			}
		}
	}

	//	---------------------------------------
	//	collision step:
	//	---------------------------------------

	for (int i=0; i<nxyzn; i++) {
		double fdiff = fprev[i] - feq[i];
		omega[i] = -fdiff/tau;
	}

	//	---------------------------------------
	//	streaming step:
	//	---------------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				for (int n=0; n<nn; n++) {
					int ndxn = (i*deli + j*delj + k*delk)*nn + n;
					int inbr = nbrIndex(i,int(ex[n]),nx);
					int jnbr = nbrIndex(j,int(ey[n]),ny);
					int nbrn = (inbr*deli + jnbr*delj + k*delk)*nn + n;
					f[nbrn]  = fprev[ndxn] + omega[ndxn];
				}
			}
		}
	}

}



// -------------------------------------------------------------------------
// Write output for the lattice-boltzmann method:
// -------------------------------------------------------------------------

void CHBasic::outputLatticeBoltzmann()
{
	int iskip = p.iskip;
	int jskip = p.jskip;
	int kskip = p.kskip;
	c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}
