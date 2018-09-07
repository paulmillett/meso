
# include "LBfluid2D.hpp"
# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// static variable initialization:
// -------------------------------------------------------------------------

int LBfluid2D::instance_count = 0;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

LBfluid2D::LBfluid2D(const CommonParams& pin) : p(pin)
{

	// ---------------------------------------
	// Unpack some of the 'params' data:
	// ---------------------------------------

	rank = p.rank;
	NX = p.NX;
	NY = p.NY;
	nx = p.nx;
	ny = p.ny;
	nxy = nx*ny;
	xOffset = p.xOff;
	nbrL = p.nbrL;
	nbrR = p.nbrR;

	// ---------------------------------------
	// Establish array dimensions:
	// ---------------------------------------

	instance_count++;
	tag = instance_count;   // instance identifier
	gx = nx + 2;            // local x-dim. + ghost nodes
	gy = ny + 2;            // local y-dim. + ghost nodes
	gxy = gx*gy;            // total lattice size
	deli = gy;              // index offset for neighbors in x-dim.
	delj = 1;               // index offset for neighbors in y-dim.

	for (int i=0; i<gxy; i++) {
		r.push_back(0.0);
		u.push_back(0.0);
		v.push_back(0.0);
		fx.push_back(0.0);
		fy.push_back(0.0);
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBfluid2D::~LBfluid2D()
{

}



// -------------------------------------------------------------------------
// Allocate f and fstream arrays (depends on stencil type):
// -------------------------------------------------------------------------

void LBfluid2D::allocateFs(const Stencil& s)
{
	int size = gxy*s.nn;
	for (int i=0; i<size; i++) {
		f.push_back(0.0);
		fstream.push_back(0.0);
	} 
}



// -------------------------------------------------------------------------
// Setters:
// -------------------------------------------------------------------------

void LBfluid2D::setTau(double val)
{
	tau = val;
}

void LBfluid2D::setRho(int i, double val)
{
	r[i] = val;
}

void LBfluid2D::setFx(int i, double val)
{
	fx[i] = val;
}

void LBfluid2D::setFy(int i, double val)
{
	fy[i] = val;
}

void LBfluid2D::setU(int i, double val)
{
	u[i] = val;
}

void LBfluid2D::setV(int i, double val)
{
	v[i] = val;
}



// -------------------------------------------------------------------------
// Getters:
// -------------------------------------------------------------------------

double LBfluid2D::getRho(int i) const
{
	return r[i];
}

double LBfluid2D::getURhoDivTau(int i) const
{
	return u[i]*r[i]/tau;
}

double LBfluid2D::getVRhoDivTau(int i) const
{
	return v[i]*r[i]/tau;
}

double LBfluid2D::getRhoDivTau(int i) const
{
	return r[i]/tau;
}



// -------------------------------------------------------------------------
// Set 'f' to 'feq'... this is ONLY done during initialization:
// -------------------------------------------------------------------------

void LBfluid2D::setFtoFeq(const Stencil& s)
{	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			int ndx = i*deli + j*delj;
			double ueq = u[ndx] + tau*fx[ndx]/r[ndx];
			double veq = v[ndx] + tau*fy[ndx]/r[ndx];
			double uv2 = ueq*ueq + veq*veq;
			for (int n=0; n<s.nn; n++) {
				int ndxn = ndx*s.nn + n;
				double evel = s.ex[n]*ueq + s.ey[n]*veq;
				double feq = 1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uv2;
				feq *= r[ndx]*s.wa[n];
				f[ndxn] = feq;
			}
		}
	}
}



// -------------------------------------------------------------------------
// Update macro arrays: 'u', 'v', 'rho':
// -------------------------------------------------------------------------

void LBfluid2D::macros(const Stencil& s, const bool ghostRho)
{
    
	// ---------------------------------------
	// update macros:
	// ---------------------------------------
	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			int ndx = i*deli + j*delj;
			double sum  = 0.0;
			double sumx = 0.0;
			double sumy = 0.0;
			for (int n=0; n<s.nn; n++) {
				int ndxn = ndx*s.nn + n;
				sum  += f[ndxn];
				sumx += f[ndxn]*s.ex[n];
				sumy += f[ndxn]*s.ey[n];
			}
			r[ndx] = sum;
			u[ndx] = sumx/r[ndx];
			v[ndx] = sumy/r[ndx];
		}
	}
	
	// ---------------------------------------
	// update rho on ghost nodes:
	// ---------------------------------------
	
	if (ghostRho) ghostNodesRho();
	
}



// -------------------------------------------------------------------------
// Streaming step:
// -------------------------------------------------------------------------

void LBfluid2D::collideStreamUpdate(const Stencil& s)
{

	// -----------------------------------
	// collision step:
	// -----------------------------------
	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			int ndx = i*deli + j*delj;
			double ueq = u[ndx] + tau*fx[ndx]/r[ndx];
			double veq = v[ndx] + tau*fy[ndx]/r[ndx];
			double uv2 = ueq*ueq + veq*veq;
			for (int n=0; n<s.nn; n++) {
				int ndxn = ndx*s.nn + n;
				double evel = s.ex[n]*ueq + s.ey[n]*veq;
				double feq = 1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uv2;
				feq *= r[ndx]*s.wa[n];
				fstream[ndxn] = f[ndxn] - (f[ndxn] - feq)/tau;
			}
		}
	}
	
	// -----------------------------------
	// exchange ghost nodes:
	// -----------------------------------

	ghostNodesStreaming(s);

	// -----------------------------------
	// streaming step:
	// -----------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			int ndx = i*deli + j*delj;
			for (int n=0; n<s.nn; n++) {
				int ndxn = (ndx)*s.nn + n;
				int inbr = i - s.exi[n];
				int jnbr = j - s.eyi[n];
				int nbrn = (inbr*deli + jnbr*delj)*s.nn + n;
				f[ndxn]  = fstream[nbrn];
			}
		}
	}
	
}



// -------------------------------------------------------------------------
// Write rho values to 'vtk' file:
// -------------------------------------------------------------------------

void LBfluid2D::writeVTKFile(std::string tagname, int tagnum,
                             int iskip, int jskip)
{

	// -----------------------------------
	//	Define the file location and name:
	// -----------------------------------

	ofstream outfile;
	std::stringstream filenamecombine;
	filenamecombine << "vtkoutput/" << tagname << "_" << tagnum << ".vtk";
	string filename = filenamecombine.str();
	outfile.open(filename.c_str(), ios::out | ios::app);

	// -----------------------------------
	//	Write the 'vtk' file header:
	// -----------------------------------

	if (rank == 0) {
		string d = "   ";
		outfile << "# vtk DataFile Version 3.1" << endl;
		outfile << "VTK file containing grid data" << endl;
		outfile << "ASCII" << endl;
		outfile << " " << endl;
		outfile << "DATASET STRUCTURED_POINTS" << endl;
		outfile << "DIMENSIONS" << d << NX/iskip << d << NY/jskip << d << 1 << endl;
		outfile << "ORIGIN " << d << 0 << d << 0 << d << 0 << endl;
		outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0 << endl;
		outfile << " " << endl;
		outfile << "POINT_DATA " << (NX/iskip)*(NY/jskip) << endl;
		outfile << "SCALARS " << tagname << " float" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
	}

	MPI::COMM_WORLD.Barrier();

	// -----------------------------------
	// Write the data:
	// NOTE: x-data increases fastest,
	//       then y-data
	// -----------------------------------

	int np = MPI::COMM_WORLD.Get_size();    // # of processors

	for (int j=1; j<ny+1; j+=jskip) {
		for (int p=0; p<np; p++) {
			if (p == rank) {
				for (int i=1; i<nx+1; i++) {
					int ig = i + xOffset;
					if (ig == 0 || ig%iskip == 0) {
						int ndx = j*delj + i*deli;
						outfile << fixed << setprecision(3) << r[ndx] << endl;
					}
				}
			}
			MPI::COMM_WORLD.Barrier();
		}
	}

	// -----------------------------------
	//	Close the file:
	// -----------------------------------

	outfile.close();

}



// -------------------------------------------------------------------------
// Update ghost nodes for rho:
// -------------------------------------------------------------------------

void LBfluid2D::ghostNodesRho()
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	int size = gy;     // size of communication
	mpiExchange(r,size,0);

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		r[i*deli + 0*delj] = r[i*deli + ny*delj];
		r[i*deli + (ny+1)*delj] = r[i*deli + 1*delj];
	}
	
}



// -------------------------------------------------------------------------
// Update ghost nodes for fstream:
// -------------------------------------------------------------------------

void LBfluid2D::ghostNodesStreaming(const Stencil& s)
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	int size = gy*s.nn;  // size of communication
	mpiExchange(fstream,size,1);

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		for (int n=0; n<s.nn; n++) {
			fstream[(i*deli + 0*delj)*s.nn + n] = fstream[(i*deli + ny*delj)*s.nn + n];
			fstream[(i*deli + (ny+1)*delj)*s.nn + n] = fstream[(i*deli + 1*delj)*s.nn + n];
		}
	}
	
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors:
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid2D::mpiExchange(std::vector<double>& a, int size, int cNum)
{
	
	MPI::Status status;
	
    // -----------------------------------
    // Send to left, Recv from right
    // -----------------------------------
	
	int ondx = 1*size;           // out index
	int indx = (nx+1)*size;	     // in index
	int stamp = tag*20 + cNum*2; // stamp is a unique int for each comm.
    
	MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrL,stamp,
	                         &a[indx],size,MPI::DOUBLE,nbrR,stamp,status);
	
	
	// -----------------------------------
	// Send to right, Recv from left
    // -----------------------------------

	ondx = nx*size;          // out index 
	indx = 0;	             // in index
	stamp += 1;              // update the stamp
	
	MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrR,stamp,
	                         &a[indx],size,MPI::DOUBLE,nbrL,stamp,status);
	
}
