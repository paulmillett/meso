
# include "MPICommunicator.hpp"


 
// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

void MPICommunicator::MPICommunicator(const CommonParams& p)
{

    nx = p.nx;
    ny = p.ny;
    deli = ;
    delj = ;

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

void MPICommunicator::~MPICommunicator()
{

}



// -------------------------------------------------------------------------
// Update ghost boundary nodes for a 2D scalar array (passed in).
// Parallelization is done by 1D domain decomposition (along x-dir.)
// -------------------------------------------------------------------------

void MPICommunicator::ghostNodesScalarArray2D(std::vector<double>& a)
{

    // -----------------------------------
    // Set boundary conditions (y-dir.)
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        a[i*deli + 0*delj] = a[ny*delj + i*deli];
        a[i*deli + (ny+1)*delj] = a[1*delj + i*deli];
    }

    // -----------------------------------
    // Set boundary conditions (x-dir.)
    // {This needs MPI communication}
    // -----------------------------------

    MPI::Status status;
    int size = deli;  // size of y-z face that will be sent

    // Send to left, Recv from right
    int stamp = tag*10;       // stamp is a unique int for each comm.
    int ondx = 1*deli;        // out index
    int indx = (nx+1)*deli;   // in  index
    MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrL,stamp,
            &a[indx],size,MPI::DOUBLE,nbrR,stamp,status);

    // Send to right, Recv from left
    stamp += 1;               // update the stamp
    ondx = nx*deli;           // out index
    indx = 0;                 // in  index
    MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrR,stamp,
            &a[indx],size,MPI::DOUBLE,nbrL,stamp,status);

}



// -------------------------------------------------------------------------
// Update ghost boundary nodes for a 2D PDF array (passed in).
// Parallelization is done by 1D domain decomposition (along x-dir.)
// -------------------------------------------------------------------------

void MPICommunicator::ghostNodesStreaming2D(std::vector<double>& a)
{

    // -----------------------------------
    // Stream top y-face to bottom y-face:
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        int ndxA = i*deli + (ny+1)*delj;
        int ndxB = i*deli + (1)*delj;
        for (int n=0; n<nn; n++) {
            if (stencil.eyi[n] > 0) {
                f[ndxB*nn + n] = f[ndxA*nn + n];
            }
        }
    }

    // -----------------------------------
    // Stream bottom y-face to top y-face:
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        int ndxA = i*deli + (0)*delj;
        int ndxB = i*deli + (ny)*delj;
        for (int n=0; n<nn; n++) {
            if (stencil.eyi[n] < 0) {
                f[ndxB*nn + n] = f[ndxA*nn + n];
            }
        }
    }

}
