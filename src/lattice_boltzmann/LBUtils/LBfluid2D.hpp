# ifndef LBFLUID2D_H
# define LBFLUID2D_H

# include "../../utils/CommonParams.h"
# include <vector>
# include <string>
# include <mpi.h>
# include "D2Q9.hpp"


class LBfluid2D {

private:

    const CommonParams& p;
    const D2Q9& stencil;
    static int instance_count;
    int tag;
    int nn;
    int NX,NY;
    int nx,ny;
    int gx,gy;
    int nxy,gxy;
    int gxyn;
    int deli,delj;
    int rank,np;
    int xOffset;
    int nbrL,nbrR;
	double tau;
	std::vector<double> u,v,r;
	std::vector<double> fx,fy;
	std::vector<double> f,feq,fstream;
	void ghostNodesStreaming();
	void mpiExchangeStreaming();
	void ghostNodesRho();
	void mpiExchangeRho();

public:

	LBfluid2D(const CommonParams&);
	~LBfluid2D();
	void macros();
    void equilDist();
	void streaming();
    void setRho(int,double);
    void setFx(int,double);
    void setFy(int,double);
	void setFeq(int,double);
	double getRho(int) const;
    void writeVTKFile(std::string,int,int,int,int);

};

# endif  // LBFLUID2D_H
