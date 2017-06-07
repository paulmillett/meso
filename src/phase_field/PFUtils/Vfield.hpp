
# ifndef VFIELD_H
# define VFIELD_H

# include "../../utils/CommonParams.h"
# include "Sfield.hpp"

class Vfield {

private:

   const CommonParams& p;
   int nx;
   int ny;
   int nz;
   int nxyz;
   Sfield ax;
   Sfield ay;
   Sfield az;

public:

	Vfield(const CommonParams&);
	~Vfield();
   void setXValues(const Sfield&);
   void setYValues(const Sfield&);
   void setZValues(const Sfield&);
   void resetVfield(std::string);
   void rescale(double);
   void addExtrapolation(double,double,double,double,double,double);
   double interpolateX(double,double,double) const;
   double interpolateY(double,double,double) const;
   double interpolateZ(double,double,double) const;
   void fft(const fftw_plan&);
   void ifft(const fftw_plan&);
	void writeVTKFile(std::string,int,int,int,int);
   Vfield& operator+=(const Vfield&);
   Vfield& operator+=(double);
   Vfield& operator-=(const Vfield&);
   Vfield& operator-=(double);
   Vfield& operator*=(const Vfield&);
   Vfield& operator*=(double);
   Vfield& operator/=(const Vfield&);
   Vfield& operator/=(double);
   Vfield& operator=(const Vfield&);
   Vfield  operator+(const Vfield&) const;
   Vfield  operator+(double) const;
   Vfield  operator-(const Vfield&) const;
   Vfield  operator-(double) const;
   Vfield  operator*(const Vfield&) const;
   Vfield  operator*(double) const;
   Vfield  operator*(fftw_complex) const;
   Vfield  operator/(const Vfield&) const;
   Vfield  operator/(double) const;
   Vfield  operator*(const Sfield&) const;
   Vfield  operator/(const Sfield&) const;
	Sfield dot(const Vfield&) const;
	Sfield div(const Vfield&) const;

};

// Some non-member methods...
const Vfield operator+(double, const Vfield&);
const Vfield operator-(double, const Vfield&);
const Vfield operator*(double, const Vfield&);
const Vfield operator*(fftw_complex, const Vfield&);
const Vfield operator*(const Sfield&, const Vfield&);

# endif  // VFIELD_H
