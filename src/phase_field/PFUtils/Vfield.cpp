
# include "Vfield.hpp"
# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

Vfield::Vfield(const CommonParams& pin) : p(pin), ax(p), ay(p), az(p)
{
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    nxyz = nx*ny*nz;
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

Vfield::~Vfield()
{

}



// -------------------------------------------------------------------------
// Setters, Getter, & Adder:
// -------------------------------------------------------------------------

void Vfield::setXValues(const Sfield& rhs)
{
    ax = rhs;
}

void Vfield::setYValues(const Sfield& rhs)
{
    ay = rhs;
}

void Vfield::setZValues(const Sfield& rhs)
{
    az = rhs;
}

Sfield* Vfield::getZValues()
{
    return &az;
}



// -------------------------------------------------------------------------
// Reset field to zero:
// -------------------------------------------------------------------------

void Vfield::resetVfield(std::string space)
{
    ax.resetSfield(space);
    ay.resetSfield(space);
    az.resetSfield(space);
}


// -------------------------------------------------------------------------
// Add a point vector to the vector field by extrpolating its value
// to neighboring grid nodes:
// -------------------------------------------------------------------------

void Vfield::addExtrapolation(double vx, double vy, double vz,
        double x, double y, double z)
{

    // determine the neighboring nodes:
    int x0 = int(floor(x)/p.dx);
    int x1 = x0 + 1;
    if (x1 >= nx) x1 = 0;
    int y0 = int(floor(y)/p.dy);
    int y1 = y0 + 1;
    if (y1 >= ny) y1 = 0;
    int z0 = int(floor(z)/p.dz);
    int z1 = z0 + 1;
    if (z1 >= nz) z1 = 0;
    double xd = x - double(x0);
    double yd = y - double(y0);
    double zd = z - double(z0);

    // --------------------------------------
    // extrapolate vx,vy,vz to those nodes:
    // --------------------------------------

    // ax.addValue(x0*nz*ny + y0*nz + z0, vx*(1.0-xd)*(1.0-yd)*(1.0-zd));
    // ay.addValue(x0*nz*ny + y0*nz + z0, vy*(1.0-xd)*(1.0-yd)*(1.0-zd));
    // az.addValue(x0*nz*ny + y0*nz + z0, vz*(1.0-xd)*(1.0-yd)*(1.0-zd));
    //
    // ax.addValue(x1*nz*ny + y0*nz + z0, vx*(xd)*(1.0-yd)*(1.0-zd));
    // ay.addValue(x1*nz*ny + y0*nz + z0, vy*(xd)*(1.0-yd)*(1.0-zd));
    // az.addValue(x1*nz*ny + y0*nz + z0, vz*(xd)*(1.0-yd)*(1.0-zd));
    //
    // ax.addValue(x0*nz*ny + y1*nz + z0, vx*(1.0-xd)*(yd)*(1.0-zd));
    // ay.addValue(x0*nz*ny + y1*nz + z0, vy*(1.0-xd)*(yd)*(1.0-zd));
    // az.addValue(x0*nz*ny + y1*nz + z0, vz*(1.0-xd)*(yd)*(1.0-zd));
    //
    // ax.addValue(x0*nz*ny + y0*nz + z1, vx*(1.0-xd)*(1.0-yd)*(zd));
    // ay.addValue(x0*nz*ny + y0*nz + z1, vy*(1.0-xd)*(1.0-yd)*(zd));
    // az.addValue(x0*nz*ny + y0*nz + z1, vz*(1.0-xd)*(1.0-yd)*(zd));
    //
    // ax.addValue(x1*nz*ny + y1*nz + z0, vx*(xd)*(yd)*(1.0-zd));
    // ay.addValue(x1*nz*ny + y1*nz + z0, vy*(xd)*(yd)*(1.0-zd));
    // az.addValue(x1*nz*ny + y1*nz + z0, vz*(xd)*(yd)*(1.0-zd));
    //
    // ax.addValue(x1*nz*ny + y0*nz + z1, vx*(xd)*(1.0-yd)*(zd));
    // ay.addValue(x1*nz*ny + y0*nz + z1, vy*(xd)*(1.0-yd)*(zd));
    // az.addValue(x1*nz*ny + y0*nz + z1, vz*(xd)*(1.0-yd)*(zd));
    //
    // ax.addValue(x0*nz*ny + y1*nz + z1, vx*(1.0-xd)*(yd)*(zd));
    // ay.addValue(x0*nz*ny + y1*nz + z1, vy*(1.0-xd)*(yd)*(zd));
    // az.addValue(x0*nz*ny + y1*nz + z1, vz*(1.0-xd)*(yd)*(zd));
    //
    // ax.addValue(x1*nz*ny + y1*nz + z1, vx*(xd)*(yd)*(zd));
    // ay.addValue(x1*nz*ny + y1*nz + z1, vy*(xd)*(yd)*(zd));
    // az.addValue(x1*nz*ny + y1*nz + z1, vz*(xd)*(yd)*(zd));

    // --------------------------------------
    // extrapolate vx,vy,vz to those nodes (evenly extrapolate):
    // --------------------------------------

    ax.addValue(x0*nz*ny + y0*nz + z0, vx*0.125);
    ay.addValue(x0*nz*ny + y0*nz + z0, vy*0.125);
    az.addValue(x0*nz*ny + y0*nz + z0, vz*0.125);

    ax.addValue(x1*nz*ny + y0*nz + z0, vx*0.125);
    ay.addValue(x1*nz*ny + y0*nz + z0, vy*0.125);
    az.addValue(x1*nz*ny + y0*nz + z0, vz*0.125);

    ax.addValue(x0*nz*ny + y1*nz + z0, vx*0.125);
    ay.addValue(x0*nz*ny + y1*nz + z0, vy*0.125);
    az.addValue(x0*nz*ny + y1*nz + z0, vz*0.125);

    ax.addValue(x0*nz*ny + y0*nz + z1, vx*0.125);
    ay.addValue(x0*nz*ny + y0*nz + z1, vy*0.125);
    az.addValue(x0*nz*ny + y0*nz + z1, vz*0.125);

    ax.addValue(x1*nz*ny + y1*nz + z0, vx*0.125);
    ay.addValue(x1*nz*ny + y1*nz + z0, vy*0.125);
    az.addValue(x1*nz*ny + y1*nz + z0, vz*0.125);

    ax.addValue(x1*nz*ny + y0*nz + z1, vx*0.125);
    ay.addValue(x1*nz*ny + y0*nz + z1, vy*0.125);
    az.addValue(x1*nz*ny + y0*nz + z1, vz*0.125);

    ax.addValue(x0*nz*ny + y1*nz + z1, vx*0.125);
    ay.addValue(x0*nz*ny + y1*nz + z1, vy*0.125);
    az.addValue(x0*nz*ny + y1*nz + z1, vz*0.125);

    ax.addValue(x1*nz*ny + y1*nz + z1, vx*0.125);
    ay.addValue(x1*nz*ny + y1*nz + z1, vy*0.125);
    az.addValue(x1*nz*ny + y1*nz + z1, vz*0.125);

}





// -------------------------------------------------------------------------
// Interpolate field value at a point between grid nodes:
// -------------------------------------------------------------------------

double Vfield::interpolateX(double x, double y, double z) const
{
    return ax.interpolate(x,y,z);
}

double Vfield::interpolateY(double x, double y, double z) const
{
    return ay.interpolate(x,y,z);
}

double Vfield::interpolateZ(double x, double y, double z) const
{
    return az.interpolate(x,y,z);
}



// -------------------------------------------------------------------------
// Rescale velocities so that maximum velocity does not exceed 'umax':
// -------------------------------------------------------------------------

void Vfield::rescale(double umax)
{

    // -----------------------------------
    //	Find largest velocity (squared):
    // -----------------------------------

    double umax2_now = 0.0;
    for (int i=0; i<nxyz; i++) {
        double ux = creal(ax.getValue(i));
        double uy = creal(ay.getValue(i));
        double uz = creal(az.getValue(i));
        double u2 = ux*ux + uy*uy + uz*uz;
        if (u2 > umax2_now) umax2_now = u2;
    }
    double umax_now_local = sqrt(umax2_now);

    // -----------------------------------
    //	Compare with all MPI processes:
    // -----------------------------------

    double umax_now = 0.0;
    MPI::COMM_WORLD.Allreduce(&umax_now_local,&umax_now,1,MPI::DOUBLE,MPI::MAX);

    cout << "Maximum velocity is " << umax_now << endl;

    // -----------------------------------
    //	Compare with 'umax' sent as input:
    // -----------------------------------

    if (umax_now > umax) {
        double rescale = umax/umax_now;

        cout << "Rescale is " << rescale << endl;

        ax = rescale*ax;
        ay = rescale*ay;
        az = rescale*az;
    }

}



// -------------------------------------------------------------------------
// FFTw transforms:
// -------------------------------------------------------------------------

void Vfield::fft(const fftw_plan& p_forward)
{
    ax.fft(p_forward);
    ay.fft(p_forward);
    az.fft(p_forward);
}

void Vfield::ifft(const fftw_plan& p_backward)
{
    ax.ifft(p_backward);
    ay.ifft(p_backward);
    az.ifft(p_backward);
}



// -------------------------------------------------------------------------
// Write array values to 'vtk' file:
// -------------------------------------------------------------------------

void Vfield::writeVTKFile(std::string tagname, int tagnum,
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

    if (p.rank == 0) {
        string d = "   ";
        outfile << "# vtk DataFile Version 3.1" << endl;
        outfile << "VTK file containing grid data" << endl;
        outfile << "ASCII" << endl;
        outfile << " " << endl;
        outfile << "DATASET STRUCTURED_POINTS" << endl;
        outfile << "DIMENSIONS" << d << p.NX/iskip << d << p.NY/jskip << d << p.NZ/kskip << endl;
        outfile << "ORIGIN " << d << 0 << d << 0 << d << 0 << endl;
        outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0*kskip << endl;
        outfile << " " << endl;
        outfile << "POINT_DATA " << (p.NX/iskip)*(p.NY/jskip)*(p.NZ/kskip) << endl;
        outfile << "VECTORS " << tagname << " float" << endl;
    }

    MPI::COMM_WORLD.Barrier();

    // -----------------------------------
    //	Write the data:
    // NOTE: x-data increases fastest,
    //       then y-data, then z-data
    // -----------------------------------

    int np = MPI::COMM_WORLD.Get_size();    // # of processors

    for (int k=0; k<p.nz; k+=kskip) {
        for (int j=0; j<p.ny; j+=jskip) {
            for (int r=0; r<p.np; r++) {
                if (r == p.rank) {
                    for (int i=0; i<p.nx; i++) {
                        int ig = i + p.xOff;
                        if (ig == 0 || ig%iskip == 0) {
                            outfile << fixed << setprecision(3)
                                << creal(ax.getValue(i*p.ny*p.nz + j*p.nz + k)) << " "
                                << creal(ay.getValue(i*p.ny*p.nz + j*p.nz + k)) << " "
                                << creal(az.getValue(i*p.ny*p.nz + j*p.nz + k)) << " "
                                << endl;
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



// -------------------------------------------------------------------------
// Compound Assignment Operators:
// -------------------------------------------------------------------------

Vfield& Vfield::operator+=(const Vfield& rhs)
{
    ax += rhs.ax;
    ay += rhs.ay;
    az += rhs.az;
    return *this;
}

Vfield& Vfield::operator+=(double val)
{
    ax += val;
    ay += val;
    az += val;
    return *this;
}

Vfield& Vfield::operator-=(const Vfield& rhs)
{
    ax -= rhs.ax;
    ay -= rhs.ay;
    az -= rhs.az;
    return *this;
}

Vfield& Vfield::operator-=(double val)
{
    ax -= val;
    ay -= val;
    az -= val;
    return *this;
}

Vfield& Vfield::operator*=(const Vfield& rhs)
{
    ax *= rhs.ax;
    ay *= rhs.ay;
    az *= rhs.az;
    return *this;
}

Vfield& Vfield::operator*=(double val)
{
    ax *= val;
    ay *= val;
    az *= val;
    return *this;
}

Vfield& Vfield::operator/=(const Vfield& rhs)
{
    ax /= rhs.ax;
    ay /= rhs.ay;
    az /= rhs.az;
    return *this;
}

Vfield& Vfield::operator/=(double val)
{
    ax /= val;
    ay /= val;
    az /= val;
    return *this;
}

Vfield& Vfield::operator=(const Vfield& rhs)
{
    ax = rhs.ax;
    ay = rhs.ay;
    az = rhs.az;
    return *this;
}



// -------------------------------------------------------------------------
// Vector Calculus:
// -------------------------------------------------------------------------

Sfield Vfield::dot(const Vfield& rhs) const
{
    return (ax*rhs.ax + ay*rhs.ay + az*rhs.az);
}

Sfield Vfield::div(const Vfield& rhs) const
{
    return ax*rhs.ax + ay*rhs.ay + az*rhs.az;
}



// -------------------------------------------------------------------------
// Binary Operators:
// -------------------------------------------------------------------------

Vfield Vfield::operator+(const Vfield& rhs) const
{
    Vfield result(p);
    result = *this;
    result += rhs;
    return result;
}

Vfield Vfield::operator+(double val) const
{
    Vfield result(p);
    result = *this;
    result += val;
    return result;
}

Vfield Vfield::operator-(const Vfield& rhs) const
{
    Vfield result(p);
    result = *this;
    result -= rhs;
    return result;
}

Vfield Vfield::operator-(double val) const
{
    Vfield result(p);
    result = *this;
    result -= val;
    return result;
}

Vfield Vfield::operator*(const Vfield& rhs) const
{
    Vfield result(p);
    result = *this;
    result *= rhs;
    return result;
}

Vfield Vfield::operator*(double val) const
{
    Vfield result(p);
    result = *this;
    result *= val;
    return result;
}

Vfield Vfield::operator*(fftw_complex val) const
{
    Vfield result(p);
    result.setXValues(ax*val);
    result.setYValues(ay*val);
    result.setZValues(az*val);
    return result;
}

Vfield Vfield::operator/(const Vfield& rhs) const
{
    Vfield result(p);
    result = *this;
    result /= rhs;
    return result;
}

Vfield Vfield::operator/(double val) const
{
    Vfield result(p);
    result = *this;
    result /= val;
    return result;
}

Vfield Vfield::operator*(const Sfield& rhs) const
{
    Vfield result(p);
    result.setXValues(ax*rhs);
    result.setYValues(ay*rhs);
    result.setZValues(az*rhs);
    return result;
}

Vfield Vfield::operator/(const Sfield& rhs) const
{
    Vfield result(p);
    result.setXValues(ax/rhs);
    result.setYValues(ay/rhs);
    result.setZValues(az/rhs);
    return result;
}



// -------------------------------------------------------------------------
// Some non-member methods:
// -------------------------------------------------------------------------

const Vfield operator+(double b, const Vfield& a)
{
    return a+b;
}

const Vfield operator-(double b, const Vfield& a)
{
    return a*(-1.0) + b;
}

const Vfield operator*(double b, const Vfield& a)
{
    return a*b;
}

const Vfield operator*(const Sfield& b, const Vfield& a)
{
    return a*b;
}

const Vfield operator*(fftw_complex b, const Vfield& a)
{
    return a*b;
}
