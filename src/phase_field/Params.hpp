
# ifndef PARAMS_H
# define PARAMS_H

# include "../utils/GetPot"



struct Params{

   GetPot input_params;
   int nx,ny,nz;
   int NX,NY,NZ;
   int rank;
   int np;
   int xOff;
   int iskip,jskip,kskip;
   double dx,dy,dz;
   double LX,LY,LZ;
   double dt;

};

# endif  // PARAMS_H
