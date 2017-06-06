
# ifndef SFIELD_H
# define SFIELD_H

# include "Params.hpp"
# include <complex.h>
# include <fftw3-mpi.h>


class Sfield {

private:

   const Params& p;
   fftw_complex* a;
   int nx,NX;
	int ny,NY;
	int nz,NZ;
   int nxyz;
   int rank,np;
   int xOffset;
   bool in_fspace;
   fftw_plan plan1;
   fftw_plan plan2;
   ptrdiff_t locsize, locnx, offx;

public:

	Sfield(const Params&);
	~Sfield();
   void setValue(int,double);
   void setValue(int,fftw_complex);
   fftw_complex getValue(int) const;
   void addValue(int,double);
   void resetSfield(std::string);
   double interpolate(double,double,double) const;
   void extrapolatePointToGrid(double,double,double,int);
   void fft(const fftw_plan&);
   void ifft(const fftw_plan&);
	void writeVTKFile(std::string,int,int,int,int);
   Sfield& operator+=(const Sfield&);
   Sfield& operator+=(double);
   Sfield& operator-=(const Sfield&);
   Sfield& operator-=(double);
   Sfield& operator*=(const Sfield&);
   Sfield& operator*=(double);
   Sfield& operator*=(fftw_complex);
   Sfield& operator/=(const Sfield&);
   Sfield& operator/=(double);
   Sfield& operator=(const Sfield&);
   Sfield  operator+(const Sfield&) const;
   Sfield  operator+(double) const;
   Sfield  operator-(const Sfield&) const;
   Sfield  operator-(double) const;
   Sfield  operator*(const Sfield&) const;
   Sfield  operator*(double) const;
   Sfield  operator*(fftw_complex) const;
   Sfield  operator/(const Sfield&) const;
   Sfield  operator/(double) const;

};

// Some non-member methods...
const Sfield operator+(double, const Sfield&);
const Sfield operator-(double, const Sfield&);
const Sfield operator*(double, const Sfield&);

# endif  // SFIELD_H
