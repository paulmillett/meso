
# ifndef FIPIAPP_H
# define FIPIAPP_H

# include "../base/MesoBase.hpp"
# include "Params.hpp"
# include "FIPISystem.hpp"

class FIPIApp : public MesoBase {

private:

	int current_step;
	Params p;
	FIPISystem* fipi_system;

public:

	FIPIApp(const GetPot&);
	~FIPIApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // FIPIAPP_H
