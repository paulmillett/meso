
# include "LBSystem.hpp"
# include <string>
# include <iomanip>
# include <sstream>
# include <fstream>
# include <iostream>
using namespace std;



// -------------------------------------------------------------------------
// Write output files:
// -------------------------------------------------------------------------

void LBSystem::writeVTKFile(string tagname, int tagnum,
                           int iskip, int jskip, int kskip)
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
   	outfile << "DIMENSIONS" << d << nxGlobal/iskip << d << ny/jskip << d << nz/kskip << endl;
   	outfile << "ORIGIN " << d << 1 << d << 1 << d << 1 << endl;
   	outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0*kskip << endl;
   	outfile << " " << endl;
   	outfile << "POINT_DATA " << (nxGlobal/iskip)*(ny/jskip)*(nz/kskip) << endl;
   	outfile << "SCALARS " << tagname << " float" << endl;
   	outfile << "LOOKUP_TABLE default" << endl;
   }

   MPI::COMM_WORLD.Barrier();

   // -----------------------------------
   //	Write the data:
   // NOTE: x-data increases fastest,
   //       then y-data, then z-data
   // -----------------------------------

   int np = MPI::COMM_WORLD.Get_size();    // # of processors

   for (int k=kstr; k<kend+1; k+=kskip) {
      for (int j=jstr; j<jend+1; j+=jskip) {
         for (int r=0; r<np; r++) {
            if (r == rank) {
               for (int i=istr; i<iend+1; i++) {
                  int ig = i + xOff;
                  if (ig == 0 || ig%iskip == 0) {
                     //outfile << fixed << setprecision(3) << rhoTot[i][j][k] << endl;
                     outfile << fixed << setprecision(3) << rho[i][j][k][0] << endl;
                  }
               }
            }
            MPI::COMM_WORLD.Barrier();
         }
      }
   }

   // -----------------------------------
   //	Close the file:
   // -----------------------------------

	outfile.close();

}
