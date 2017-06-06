
# ifndef FFTARRAYS_H
# define FFTARRAYS_H

# include "Params.hpp"
# include <string>
# include <complex.h>
# include <fftw3-mpi.h>
using namespace std;

class FFTArrays {

private:

	int nxGlobal;
	int xOffset;
	int nx;
	int ny;
	int nz;
	int rank;
	int sizeR;
	int sizeC;
	double dx;
	double dy;
	double dz;
	double dt;
	double* c;
	double* mu;
	double* k1;
	double* k2;
	double* k4;
	double* mob;
	double* gradmu;
	double* Mgradmu;
	fftw_complex* cFFT;
	fftw_complex* muFFT;
	fftw_complex* gradmuFFT;
	fftw_complex* MgradmuFFT;
	fftw_complex i1;
	fftw_plan plan1;  // forward transform plan for 'c'
	fftw_plan plan2;  // forward transform plan for 'mu'
	fftw_plan plan3;  // backward transform plan for 'c'
	fftw_plan plan4;  // backward transform plan for 'gradmu'
	fftw_plan plan5;  // forward transform plan for 'Mgradmu'
	ptrdiff_t locsize, locnx, offx;

public:

	FFTArrays(const Params&);
	~FFTArrays();
	double getConc(int) const;
	void setConc(int,double);
	void addtoConc(int,double);
	void setMu(int,double);
	void addtoMu(int,double);
	void setMob(int,double);
	double getMob(int) const;
	int getRealSize();
	int getNxLocal();
	int getXoffset();
	void forwardFFT();
	void backwardFFT();
	void updateConcFFT(bool);
	void writeVTKFile(std::string,int,int,int,int);

};

# endif  // FFTARRAYS_H
