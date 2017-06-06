
# ifndef BDAPP_H
# define BDAPP_H

# include <vector>
# include "../base/MesoBase.hpp"
# include "Particles.hpp"
# include "Params.hpp"

class BDApp : public MesoBase {

private:

	int current_step;
	Params p;
	Particles* bd_particles;

public:

	BDApp(const GetPot&);
	~BDApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // BDAPP_H
