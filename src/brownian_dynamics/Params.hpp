
# ifndef PARAMS_H
# define PARAMS_H

# include "../utils/GetPot"



struct Params{

   GetPot input_params;
   int N;
   int np,rank;
   double Lx,Ly,Lz;
   double dt;
   
};

# endif  // PARAMS_H
