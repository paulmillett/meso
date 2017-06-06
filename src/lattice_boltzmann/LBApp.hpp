
# ifndef LBAPP_H
# define LBAPP_H

# include <vector>
# include "../base/MesoBase.hpp"
# include "LBSystem.hpp"
# include "mcmp.hpp"

class LBApp : public MesoBase {

private:

	int current_step;
	LBSystem* lb_system;
	mcmp* mcmp_object;

public:

	LBApp(const GetPot&);
	~LBApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // LBAPP_H
