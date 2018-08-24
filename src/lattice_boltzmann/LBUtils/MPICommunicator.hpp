# ifndef MPICOMMUNICATOR_H
# define MPICOMMUNICATOR_H

# include "../../utils/CommonParams.h"
# include <vector>
# include <string>
# include <mpi.h>


class MPICommunicator {

private:

    int nx,ny;
    int deli,delj;
    int rank,np;
    int xOffset;
    int nbrL,nbrR;

public:

	MPICommunicator(const CommonParams&);
	~MPICommunicator();
    void ghostNodesScalarArray2D(std::vector<double>&);
    void ghostNodesStreaming2D(std::vector<double>&);

};
 
# endif  // MPICOMMUNICATOR_H
