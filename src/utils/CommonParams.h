
# ifndef COMMONPARAMS_H
# define COMMONPARAMS_H

struct CommonParams{

   int nx,ny,nz;
   int NX,NY,NZ;
   int rank;
   int np;
   int nbrL;
   int nbrR;
   int xOff;
   int iskip;
   int jskip;
   int kskip;
   int nstep;
   int numOutputs;
   double dx,dy,dz;
   double LX,LY,LZ;
   double dt;

};

# endif  // COMMONPARAMS_H
