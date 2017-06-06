

# ifndef CUBICLATTICE_H
# define CUBICLATTICE_H

# include "PInitCond_Interface.hpp"

/*
   This is an implementation of the 'PInitCond' interface class.
   It initializes a set of particles on a 3D cubic lattice.
*/

class CubicLattice : public PInitCond {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int N;
   double Lx, Ly, Lz;
   double vscl;
   double pradii;
   vector<double>& r;
   vector<double>& v;
   vector<double>& rad;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

   // -----------------------------
   // Constructor...
   // -----------------------------

   CubicLattice(const Params& p, vector<double>& rin,
                vector<double>& vin, vector<double>& radin) :
                r(rin), v(vin), rad(radin)
   {
      N = p.input_params("BDApp/N",1);
      Lx = p.input_params("BDApp/Lx",5.0);
      Ly = p.input_params("BDApp/Ly",5.0);
      Lz = p.input_params("BDApp/Lz",5.0);
      vscl = p.input_params("BDApp/initial_condition/vscl",0.0);
      pradii = p.input_params("BDApp/initial_condition/pradii",1.0);
   }


   // -----------------------------
   // Destructor...
   // -----------------------------

   ~CubicLattice()
   {
   }

   // -----------------------------
   // Function to calculate i.c.:
   // -----------------------------

   void icFunc()
   {

      // initialize particles in a cubic array:
      int nx = 10;
      int ny = nx;
      int nz = nx;
      double dr[3] = {Lx/nx, Ly/ny, Lz/nz};
      double ic[3] = {0.0, 0.0, 0.0};
      int parti = 0;
      for (int ix=0; ix<nx; ix++) {
         ic[0]=ix;
         for (int iy=0; iy<ny; iy++) {
            ic[1]=iy;
            for (int iz=0; iz<nz; iz++) {
               ic[2]=iz;
               for (int k=0; k<3; k++) {
                  r[parti*3+k] = dr[k]*(double(ic[k])+0.5);
               }
               parti++;
            }
         }
      }

      // initialize particle radii:
      for (int i=0; i<N; i++) {
         rad[i] = pradii;
      }

      // initialize particle velocities:
      double vsum[3] = {0.0, 0.0, 0.0};
      for (int i=0; i<N; i++) {
         for (int k=0; k<3; k++) {
            double r = (double)rand()/RAND_MAX;
            v[i*3+k] = vscl*(r - 0.5);
            vsum[k] += v[i*3+k];
         }
      }

      // zero the total momentum:
      for (int i=0; i<N; i++) {
         for (int k=0; k<3; k++) {
            v[i*3+k] -= vsum[k]/N;
         }
      }

   }

};

# endif  // CUBICLATTICE_H
