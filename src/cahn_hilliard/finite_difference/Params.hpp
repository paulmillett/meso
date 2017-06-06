
# ifndef PARAMS_H
# define PARAMS_H

# include "../../utils/GetPot"



struct Params{

   GetPot input_params;
   int nx,NX;
   int ny,NY;
   int nz,NZ;
   int nc;
   int rank,nbrL,nbrR;
   int np;
   int xOff;
   double dx,dy,dz;
   double dt;

};

# endif  // PARAMS_H
