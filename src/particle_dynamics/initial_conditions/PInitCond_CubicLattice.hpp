

# ifndef CUBICLATTICE_H
# define CUBICLATTICE_H

# include "PInitCond_Interface.hpp"

/*
   This is an implementation of the 'PInitCond' interface class.
   It initializes a set of particles on a 3D cubic lattice. It
   is assume that the number of particles has an integral cube
   root i.e. (N)^(1/3) is an integer.
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

   CubicLattice(const GetPot& p, vector<double>& rin,
                vector<double>& vin, vector<double>& radin) :
                r(rin), v(vin), rad(radin)
   {
      N = p("BDApp/N",1);
      Lx = p("BDApp/Lx",5.0);
      Ly = p("BDApp/Ly",5.0);
      Lz = p("BDApp/Lz",5.0);
      vscl = p("BDApp/initial_condition/vscl",0.0);
      pradii = p("BDApp/initial_condition/pradii",1.0);
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
      int nx = pow((double)N+1,1.0/3.0); //number of particles must have integral cube root
      int ny = nx;
      int nz = nx;
      int part_index;
      double a0 = (double)Lx/(double)nx;
      double xpos,ypos,zpos;
      double r1,r2,r3;
      for (int k=0; k<nx; k++)
      {
          for (int j=0; j<nx; j++)
          {
              for (int i=0; i<nx; i++)
              {
                  r1 = (double)rand()/RAND_MAX-0.5;
                  r3 = (double)rand()/RAND_MAX-0.5;
                  r2 = (double)rand()/RAND_MAX-0.5;
                  part_index = k*nx*nx+j*nx+i;
                  xpos = a0/2.0 + (double)i*a0;
                  ypos = a0/2.0 + (double)j*a0;
                  zpos = a0/2.0 + (double)k*a0;
                  r[3*part_index+0] = zpos + 0.3*r1*a0;
                  r[3*part_index+1] = ypos + 0.3*r2*a0;
                  r[3*part_index+2] = xpos + 0.3*r3*a0;

              }
          }
      }
//      double dr[3] = {Lx/nx, Ly/ny, Lz/nz};
//      double ic[3] = {0.0, 0.0, 0.0};
//      int parti = 0;
//      for (int ix=0; ix<nx; ix++) {
//         ic[0]=ix;
//         for (int iy=0; iy<ny; iy++) {
//            ic[1]=iy;
//            for (int iz=0; iz<nz; iz++) {
//               ic[2]=iz;
//               for (int k=0; k<3; k++) {
//                  r[parti*3+k] = dr[k]*(double(ic[k])+0.5);
//               }
//               parti++;
//            }
//         }
//      }

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
