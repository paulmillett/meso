
# ifndef CHFFTAPP_H
# define CHFFTAPP_H

# include <vector>
# include "../../base/MesoBase.hpp"
# include "CHFFTSystem.hpp"
# include "Params.hpp"

class CHFFTApp : public MesoBase {

private:

	int current_step;
	Params p;
	CHFFTSystem* ch_system;

public:

	CHFFTApp(const GetPot&);
	~CHFFTApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // CHFFTAPP_H
