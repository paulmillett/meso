
# ifndef PARAMS_H
# define PARAMS_H

# include "../utils/GetPot"



struct Params{

   int N;
   int nx,NX;
   int ny,NY;
   int nz,NZ;
   int rank;
   int np;
   int xOff;
   int dtRatio;
   int iskip,jskip,kskip;
   double dx,dy,dz;
   double Lx,Ly,Lz;
   double dt;
   double visc;
   double pmob;
   double rcut;
   double A;
   double scl;
   double Khertz;
   double co;
   double w,kap;

};

# endif  // PARAMS_H
