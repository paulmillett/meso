
# include "LBfluid2D.hpp"
# include "GhostBoundaryNodes.hpp"
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

LBfluid2D::LBfluid2D(const CommonParams& pin,
                     const D2Q9 sin) : p(pin), stencil(sin)
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
    nn = stencil.nn;
    xOffset = p.xOff;
    nbrL = p.nbrL;
    nbrR = p.nbrR;

    // ---------------------------------------
    // Establish array dimensions:
    // ---------------------------------------

    instance_count++;
    tag = instance_count;   // array identifier
    gx = nx + 2;            // local x-dim. + ghost nodes
    gy = ny + 2;            // local y-dim. + ghost nodes
    gxy = gx*gy;            // total lattice size
    gxyn = gxy*nn;          // total size of vector arrays
    deli = gy;              // index offset for neighbors in x-dim.
    delj = 1;               // index offset for neighbors in y-dim.

    for (int i=0; i<gxy; i++) {
        r.push_back(0.0);
        u.push_back(0.0);
        v.push_back(0.0);
        fx.push_back(0.0);
        fy.push_back(0.0);
        for n=0; n<nn; n++) {
            f.push_back(0.0);
            feq.push_back(0.0);
            fstream.push_back(0.0);
        }
    }

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBfluid2D::~LBfluid2D()
{

}



// -------------------------------------------------------------------------
// Setters and Getters:
// -------------------------------------------------------------------------

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

void LBfluid2D::setFeq(int i, double val)
{
	feq[i] = val;
}

double LBfluid2D::getRho(int i) const
{
    return r[i];
}



// -------------------------------------------------------------------------
// Update macro arrays: 'u', 'v', 'rho':
// -------------------------------------------------------------------------

void LBfluid2D::macros()
{
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            int ndx = i*deli + j*delj;
            double sum  = 0.0;
            double sumx = 0.0;
            double sumy = 0.0;
            for (int n=0; n<nn; n++) {
                int ndxn = ndx*nn + n;
                sum  += f[ndxn];
                sumx += f[ndxn]*stencil.ex[n];
                sumy += f[ndxn]*stencil.ey[n];
            }
            r[ndx] = sum;
            u[ndx] = sumx/r[ndx];
            v[ndx] = sumy/r[ndx];
        }
    }
    ghostNodesRho();
}



// -------------------------------------------------------------------------
// Update equilibrium particle distributions:
// -------------------------------------------------------------------------

void LBfluid2D::equilDist()
{
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            int ndx = i*deli + j*delj;
            double u_F = u[ndx] + tau*fx[ndx]/r[ndx];
            double v_F = v[ndx] + tau*fy[ndx]/r[ndx];
            double uv2 = u_F*u_F + v_F*v_F;
            for (int n=0; n<nn; n++) {
                int ndxn = ndx*nn + n;
                double evel = stencil.ex[n]*u_F + stencil.ey[n]*v_F;
                feq[ndxn] = 1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uv2;
                feq[ndxn] *= r[ndx]*stencil.wa[n];
            }
        }
    }
}



// -------------------------------------------------------------------------
// Streaming step:
// -------------------------------------------------------------------------

void LBfluid2D::streaming()
{

    // -----------------------------------
    // collision step:
    // -----------------------------------

    for (int i=0; i<gxyn; i++) {
        double fdiff = f[i] - feq[i];
        fstream[i] = f[i] - fdiff/tau;
    }

    // -----------------------------------
    // exchange ghost nodes:
    // -----------------------------------

    ghostNodesStreaming();

    // -----------------------------------
    // streaming step:
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
			int ndx = i*deli + j*delj;
            for (int n=0; n<nn; n++) {
                int ndxn = (ndx)*nn + n;
                int inbr = i - stencil.exi[n];
                int jnbr = j - stencil.eyi[n];
                int nbrn = (inbr*deli + jnbr*delj)*nn + n;
                f[ndxn]  = fstream[nbrn];
            }
        }
    }
	
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors (fstream):
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid2D::ghostNodesStreaming()
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	mpiExchangeStreaming();

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		for (int n=0; n<nn; n++) {
			fstream[(i*deli + 0*delj)*nn + n] = fstream[(i*deli + ny*delj)*nn + n];
			fstream[(i*deli + (ny+1)*delj)*nn + n] = fstream[(i*deli + 1*delj)*nn + n];
		}
	}
	
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors (fstream):
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid2D::mpiExchangeStreaming()
{
	
	MPI::Status status;
	int size = gy*nn;
	
    // -----------------------------------
    // Send to left, Recv from right
    // -----------------------------------
		
	int ondx = 1*size;       // out index
	int indx = (nx+1)*size;	 // in index
	int stamp = tag*10;      // stamp is a unique int for each comm.
    
	MPI::COMM_WORLD.Sendrecv(&fstream[ondx],size,MPI::DOUBLE,nbrL,stamp,
	                         &fstream[indx],size,MPI::DOUBLE,nbrR,stamp,status);
	
    // -----------------------------------
    // Send to right, Recv from left
    // -----------------------------------
		
	ondx = nx*size;          // out index 
	indx = 0;	             // in index
	stamp += 1;              // update the stamp
	
	MPI::COMM_WORLD.Sendrecv(&fstream[ondx],size,MPI::DOUBLE,nbrR,stamp,
 	                         &fstream[indx],size,MPI::DOUBLE,nbrL,stamp,status);

}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors (rho):
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid2D::ghostNodesRho()
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	mpiExchangeRho();

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		r[i*deli + 0*delj] = r[i*deli + ny*delj];
		r[i*deli + (ny+1)*delj] = r[i*deli + 1*delj];
	}
	
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors (fstream):
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid2D::mpiExchangeRho()
{
	
	MPI::Status status;
	int size = gy;
	
    // -----------------------------------
    // Send to left, Recv from right
    // -----------------------------------
		
	int ondx = 1*size;       // out index
	int indx = (nx+1)*size;	 // in index
	int stamp = tag*10;      // stamp is a unique int for each comm.
    
	MPI::COMM_WORLD.Sendrecv(&r[ondx],size,MPI::DOUBLE,nbrL,stamp,
	                         &r[indx],size,MPI::DOUBLE,nbrR,stamp,status);
	
    // -----------------------------------
    // Send to right, Recv from left
    // -----------------------------------
		
	ondx = nx*size;          // out index 
	indx = 0;	             // in index
	stamp += 1;              // update the stamp
	
	MPI::COMM_WORLD.Sendrecv(&r[ondx],size,MPI::DOUBLE,nbrR,stamp,
 	                         &r[indx],size,MPI::DOUBLE,nbrL,stamp,status);

}



// -------------------------------------------------------------------------
// Write rho values to 'vtk' file:
// -------------------------------------------------------------------------

void SfieldFD::writeVTKFile(std::string tagname, int tagnum,
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
        outfile << "DIMENSIONS" << d << NX/iskip << d << NY/jskip << endl;
        outfile << "ORIGIN " << d << 0 << d << 0 << endl;
        outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << endl;
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
        for (int r=0; r<np; r++) {
            if (r == rank) {
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
